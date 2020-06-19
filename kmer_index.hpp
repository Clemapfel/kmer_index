//
// Created by Clemens Cords on 2/6/20.
// Copyright (c) 2020 Clemens Cords. All rights reserved.
//

#pragma once

#include <cmath>
#include <seqan3/alphabet/concept.hpp>
#include <seqan3/alphabet/hash.hpp>
#include <seqan3/range/views/kmer_hash.hpp>
#include <seqan3/core/type_traits/range.hpp>
#include <seqan3/core/debug_stream.hpp>

#include <type_traits>
#include <cstdint>
#include <cassert>
#include <cmath>
#include <algorithm>
#include <limits>
#include <unordered_map>

#include <robin_hood.h>
#include <fast_pow.hpp>

#include <kmer_index_result.hpp>
#include <thread_pool.hpp>
#include <compressed_bitset.hpp>


namespace kmer
{
    namespace detail
    {
        // represents a kmer index for a single k
        template<seqan3::alphabet alphabet_t, typename position_t, size_t k>
        class kmer_index_element
        {
            // k has to be != 0 and small enough to be representable with size_t as hash
            static_assert(k > 0 and k < 64 / log2(seqan3::alphabet_size<alphabet_t>), "please specify a valid k");

            private:
                // typedefs for readability
                using result_t = kmer_index_result<position_t>;
                constexpr static size_t _sigma = seqan3::alphabet_size<alphabet_t>;

                // data
                robin_hood::unordered_map<size_t, std::vector<position_t>> _data;

                // needed for subk search edge case
                std::vector<alphabet_t> _last_kmer;
                std::vector<std::vector<position_t>> _last_kmer_refs;

                /*
                // hash a query of arbitrary length
                template<typename iterator_t>
                static size_t hash_any(iterator_t query_it, size_t size)
                {
                    auto it = query_it;

                    size_t hash = 0;
                    for (size_t i = 0; i < size; ++i)
                        hash += (seqan3::to_rank(*it++) * detail::fast_pow(_sigma, size - i - 1));

                    return hash;
                }
                */

                // hash a query of length k, optimized by compiler unwrapping fold expression bc k is constexpr
                // aux_aux redundant but forces proper order of evaluation during unwrap
                template<typename iterator_t>
                size_t hash_aux_aux(iterator_t query_it, size_t i) const
                {
                    return seqan3::to_rank(*query_it) * kmer::detail::fast_pow(_sigma, k - i - 1);
                }

                template<typename iterator_t, size_t... is>
                size_t hash_aux(iterator_t query_it, std::index_sequence<is...> sequence) const
                {
                    return (... + hash_aux_aux(query_it++, is));
                }

                template<typename iterator_t>
                size_t hash(iterator_t query_it) const
                {
                    auto it = query_it;
                    return hash_aux(it, std::make_index_sequence<k>());
                }

                // get vector of positions from map
                const std::vector<position_t>* at(size_t hash) const
                {
                    auto it = _data.find(hash);

                    if (it != _data.end())
                        return &it->second;
                    else
                        return nullptr;
                }

                // check last kmer (needed for edge case where subk query is at the very end of the text)
                template<typename iterator_t>
                void check_last_kmer(iterator_t subk_begin, size_t size,
                                     std::vector<const std::vector<position_t>*>& to_fill) const
                {
                    for (size_t i = 1; i < k - size + 1; ++i)   // first char is covered by last kmer being in _data
                    {
                        bool equal = true;
                        auto it = subk_begin;

                        for (size_t j = i; j < i + size; ++j)
                        {
                            if (_last_kmer.at(j) != *it)
                            {
                                equal = false;
                                break;
                            }
                            it++;
                        }

                        if (equal)
                            to_fill.push_back(&_last_kmer_refs.at(i));
                            // returns pointer to vector with correct position at text.size() - k + i
                    }
                }

                // generate hashs for all kmer that have the given prefix and lookup those hashes
                // used for searching queries of lenght < k, a kmer with prefix p == query has the
                // same position as the query itself
                template<typename iterator_t>
                std::vector<const std::vector<position_t>*>
                get_position_for_all_kmer_with_prefix(iterator_t prefix_begin, size_t size) const
                {
                    if (fast_pow(_sigma, k-size) >= 5e5)
                    {
                        std::cerr << "[WARNING] lookup of query with size " << size << " with kmer for k = " << k
                                  << " may take a long time\n";

                        //if (fast_pow(_sigma, k-size) >= 2e7)
                            //throw std::invalid_argument("invalid query, please choose a different k");
                    }

                    auto it = prefix_begin;
                    size_t prefix_hash = 0;
                    for (size_t i = 0; i < size; ++i)
                    {
                        prefix_hash += seqan3::to_rank(*it++) * fast_pow(_sigma, k - i - 1);
                    }

                    //c.f. addendum below
                    size_t lower_bound = 0 + prefix_hash;
                    size_t n_hashes =  kmer::detail::fast_pow(_sigma, k - size);
                    size_t upper_bound = lower_bound + n_hashes;
                    size_t step_size = 1;

                    std::vector<const std::vector<position_t>*> output;

                    // bc stepsize is 1, just add number of possible hashes to lower_bound;
                    for (size_t hash = lower_bound;
                         hash < upper_bound; hash += step_size)
                    {
                        const auto* pos = at(hash);
                        if (pos != nullptr)
                            output.push_back(pos);
                    }

                    check_last_kmer(prefix_begin, size, output);
                    return output;

                    // Addendum:
                    // to calculat set H of all hashes with given prefix p
                    // i)  hash of any kmer q_0, q_1, ... q_k = sum {i=0 to k} rank(q_i) * sigma^(k - i - 1)
                    //
                    // ii) for given prefix p_0, ..., p_m, the first part of the sum is known:
                    //     h_p = sum {i=0 to m} rank(p_i) * sigma^(k - i - 1)
                    //
                    // iii) for arbitrary suffix s_0, ..., s_(k-m-1) we observe:
                    //      min(H) = h_p    by setting all s_0 to a char that has rank 0;
                    //
                    //  iv) max(H) = h_p + (hash of suffix with all s_i so that rank(s_i) = maximum = sigma-1)
                    //             = h_p + sum {i=k-m to i=k} (sigma-1) * sigma^(k-i-1)
                    //             = h_p + sigma^m - 1/sigma
                    //
                    //   v) for two hashes from H h_1 and h_2 so that h_1 < h_2 it is true that
                    //      h_1 - h_2 >= 1 as the minimum step to get from h_1 to h_2 is to set the
                    //      last char of suffix of h_2 to be one rank higher as the last char in h_1.
                    //      The last char contributes sigma^(k-k-1-1) = sigma^0 = 1 to the hash sum
                    //
                    //   Thus to generate all hashes in H we calculate the lower bound min(H) = h_p,
                    //   the upper bound max(H) = h_p + sigma^m - 1/sigma and go stepwis se by 1 from minimum
                    //   maximum. This is to the authors knowledge the fastest way to generate all hashes in H.
                    //   For the implementation max(H) = h_p + sigma^(k-m) is used, since min(H) = h_p and
                    //   #H = sigma^(k-m) and every hash is 1 higher than the previous hash
                    //   max(H) = h_p + #H = sigma^(k-m) which is one less operation than the above formula
                }

            // all functions are at least protected bc the user should only use index_element functionality through kmer_index
            public:

                // search any query, bevahior (and thus runtime) dependend on query length
                virtual result_t search(std::vector<alphabet_t>& query) const
                {
                    assert(query.size() > 0);

                    // for query.size() = m == k, just do a single lookup in the map
                    if (query.size() == k)
                    {
                        const auto* pos = at(hash(query.begin()));
                        if (pos)
                            return result_t(pos);
                        else
                            return result_t();
                    }
                    // for m > k split query into m / k (+1) parts and look up each one individually
                    else if (query.size() > k)
                    {
                        size_t rest_n = query.size() % k;

                        // get positions for nk parts
                        std::vector<const std::vector<position_t>*> nk_positions;

                        int last_hash = -1;

                        for (size_t i = 0; i < query.size() - rest_n; i += k)
                        {
                            // tiny optimization: skip map query if part has the same hash
                            size_t h = hash(query.begin() + i);
                            const auto* pos = (h != last_hash ? at(h) : nk_positions.back());

                            if (pos)
                                nk_positions.push_back(pos);
                            else
                                return result_t();

                            last_hash = h;
                        }

                        // precompute rest positions
                        auto usable = detail::compressed_bitset(nk_positions.back()->size(), true);

                        if (rest_n > 0)
                        {
                            auto rest_results = get_position_for_all_kmer_with_prefix(query.end() - rest_n, rest_n);

                            size_t i = 0;
                            for (auto pos : *nk_positions.back())
                            {
                                bool can_be_used = false;
                                for (const auto* res_vec : rest_results)
                                {
                                    if (std::binary_search(res_vec->begin(), res_vec->end(), pos + k))
                                    {
                                        can_be_used = true;
                                        break;
                                    }
                                }

                                if (can_be_used)
                                    usable.set_1(i);
                                else
                                    usable.set_0(i);

                                ++i;
                            }
                        }

                        // query.size % k != 0 and query.size < 2*k
                        if (nk_positions.size() == 1)
                        {
                            result_t output(nk_positions.at(0), false);
                            for (size_t i = 0; i < nk_positions.at(0)->size(); ++i)
                                if (usable.at(i))
                                    output.should_use(i);

                            return output;
                        }

                        // query.size % k == 0 does not need to check rest
                        else if (rest_n == 0)
                        {
                            result_t output(nk_positions.at(0), true);

                            for (size_t start_pos_i = 0; start_pos_i < nk_positions.front()->size(); ++start_pos_i)
                            {
                                size_t previous_pos = nk_positions.front()->at(start_pos_i);
                                bool should_use = true;

                                for (size_t next_pos_i = 1; next_pos_i < nk_positions.size(); ++next_pos_i)
                                {
                                    auto* current = nk_positions.at(next_pos_i);

                                    // check if next part has a position that is equal to the current + k
                                    if (not std::binary_search(current->begin(), current->end(), previous_pos + k))
                                    {
                                        output.should_not_use(start_pos_i);
                                        should_use = false;
                                        break;
                                    }
                                    else
                                        previous_pos += k;
                                }

                                if (should_use)
                                    output.should_use(start_pos_i);
                            }

                            return output;
                        }

                        // query.size() % k != 0 and query.size() > 2*k
                        else
                        {
                            result_t output(nk_positions.at(0), false);

                            for (size_t start_pos_i = 0; start_pos_i < nk_positions.front()->size(); ++start_pos_i)
                            {
                                size_t previous_pos = nk_positions.front()->at(start_pos_i);

                                if (nk_positions.size() > 1)
                                {
                                    bool interrupted = false;
                                    for (size_t next_pos_i = 1; next_pos_i < nk_positions.size(); ++next_pos_i)
                                    {
                                        const auto* current = nk_positions.back();
                                        auto it = std::lower_bound(current->begin(), current->end(), previous_pos += k);

                                        if (*it != previous_pos)
                                        {
                                            interrupted = true;
                                            break;
                                        }

                                        // for last k part, also check precomputed rest pos
                                        if (next_pos_i == nk_positions.size() - 1)
                                        {
                                            if (not usable.at(it - current->begin()))
                                                interrupted = true;

                                            break;
                                        }
                                    }

                                    if (interrupted)
                                        output.should_not_use(start_pos_i);
                                }
                            }

                            return output;
                        }
                    }

                    // for query.size() < k use prefix search directly
                    else
                    {
                        return result_t(get_position_for_all_kmer_with_prefix(query.begin(), query.size()), true);
                    }
                }

                // CTORS
                kmer_index_element() = default;

                template<std::ranges::range text_t>
                kmer_index_element(text_t& text)
                        :  _last_kmer(text.end() - k, text.end()),
                           _last_kmer_refs([this]() -> std::vector<std::vector<position_t>> {
                               std::vector<std::vector<position_t>> output;

                               for (position_t i = 0; i < _last_kmer.size(); ++i)
                                   output.emplace_back(std::vector<position_t>{i});

                               return output;
                           }())
                {
                    create(text);
                }

                // create is a seperate function so kmer_index below can run multiple creates in paralell
                template<std::ranges::range text_t>
                void create(text_t& text)
                {
                    auto hashes = text | seqan3::views::kmer_hash(seqan3::shape{seqan3::ungapped{k}});

                    size_t i = 0;
                    for (auto h : hashes)
                    {
                        if (_data.find(h) == _data.end())
                            _data.emplace(h, std::vector<position_t>());

                        _data[h].push_back(i);
                        ++i;
                    }

                    assert(i + k -1 < std::numeric_limits<position_t>::max() && "your text is too large for this configuration. Please specify position_t = uint64_t manually");

                    // setup last kmer
                    position_t text_size = position_t(i) + k - 1;

                    _last_kmer = std::vector<alphabet_t>(text.end() - k, text.end());

                    _last_kmer_refs = std::vector<std::vector<position_t>>();

                    for (position_t j = 0; j < _last_kmer.size(); ++j)
                        _last_kmer_refs.emplace_back(std::vector<position_t>{j + text_size - position_t(k)});
                }
        };
    } // end of namespace detail

    template<seqan3::alphabet alphabet_t, typename position_t, size_t... ks>
    class kmer_index
        : public detail::kmer_index_element<alphabet_t, position_t, ks>...
    {
        private:
            // typedefs for readability
            template<size_t k>
            using index_element_t = detail::kmer_index_element<alphabet_t, position_t, k>;
            using index_t = kmer_index<alphabet_t, position_t, ks...>;
            using result_t = detail::kmer_index_result<position_t>;

            // preallocate list of all ks as it is used often
            inline static std::vector<size_t> _all_ks = std::vector<size_t>{ks...};

            // search for a specific k index_element has to be able to be called thus...
            template<size_t k>
            result_t call_search(std::vector<alphabet_t>& query) const
            {
                return static_cast<const index_element_t<k>*>(this)->index_element_t<k>::search(query);
            }

            typedef result_t(kmer_index<alphabet_t, position_t, ks...>::*element_search_fn)(
                    std::vector<alphabet_t>&) const;

            // ... setup array that holds the corresponding search call so it can be a simple lookup.
            // This allows for RVO as during a fold expression the result would have to be handed from the
            // index_element to this index and it's move constructor would have to be called at least once
            const std::array<element_search_fn, sizeof...(ks)> _search_fns = {
                    (&kmer_index<alphabet_t, position_t, ks...>::call_search<ks>)...};

            // array that holds the index of _search_fns that corresponds to a specific k
            // further modified in ctor
            std::array<int, std::max({ks...}) + 1> _k_to_search_fns_i;

            // precalculate what ks to use based on query length (2000 is maximum common read length)
            inline static constexpr size_t _query_size_range = 2000;
            std::array<size_t, _query_size_range> _optimal_k;

            size_t choose_best_k_for_query_size(size_t query_size) const
            {
                // pick best k to search with
                size_t optimal_k = _all_ks.front(); // _all_ks need to be sorted highest to lowest

                if (_all_ks.size() > 1)
                {
                    for (size_t k : _all_ks)
                    {
                        if (query_size % k == 0)
                        {
                            optimal_k = k;
                            break;
                        }
                        // the next highest multiple of k > query_size is as close to query_size as possible
                        else if ((std::ceil(query_size / float(k))*k - query_size) < (std::ceil(query_size / float(optimal_k))*optimal_k - query_size))
                            optimal_k = k;
                    }
                }

                return optimal_k;
            }

        public:
            // ctor
            template<std::ranges::range text_t>
            kmer_index(text_t& text, size_t n_threads = std::thread::hardware_concurrency())
                    : index_element_t<ks>()...
            {
                // catch hardware concurrency failing
                if (n_threads == 0)
                    n_threads = 1;

                // use multiple threads to build index elements at the same time
                auto pool = detail::thread_pool{n_threads};

                std::vector<std::future<void>> futures;
                (futures.emplace_back(
                        pool.execute(&index_element_t<ks>::template create<text_t>,
                                     static_cast<index_element_t<ks>*>(this),
                                     std::ref(text))), ...);
                // wait to finish
                for (auto& f : futures)
                    f.get();

                // setup k_to_search_fn_i
                std::fill(_k_to_search_fns_i.begin(), _k_to_search_fns_i.end(), -1);
                size_t i = 0;
                for (size_t k : _all_ks)
                {
                    _k_to_search_fns_i[k] = i;
                    ++i;
                }

                // sort all_ks so bigger ks can be prioritized in search
                std::sort(_all_ks.begin(), _all_ks.end(), [](size_t a, size_t b) -> bool {return a > b;});

                for (size_t query_size = 0; query_size < _optimal_k.size(); ++query_size)
                    _optimal_k[query_size] = choose_best_k_for_query_size(query_size);
            }

            // search single query with index_element<k> directly
            template<size_t k>
            result_t search(std::vector<alphabet_t>& query) const
            {
                return index_element_t<k>::search(query);
            }

            // search single query, index picks optimal search scheme.
            // no overhead thanks to RVO
            result_t search(std::vector<alphabet_t>& query) const
            {
                // call corresponding index_element search
                // optimal k for common queries (read length <= 2000) pre-computed
                if (query.size() < _query_size_range)
                    return (this->*_search_fns[_k_to_search_fns_i.at(_optimal_k.at(query.size()))])(query);
                else
                    return (this->*_search_fns[_k_to_search_fns_i.at(choose_best_k_for_query_size(query.size()))])(query);
            }

            // search multiple queries in paralell
            template<std::ranges::forward_range queries_t>
            std::vector<result_t>
            search(queries_t&& queries, size_t n_threads = std::thread::hardware_concurrency()) const
            {
                if (n_threads == 0)
                    n_threads = 1;

                auto pool = detail::thread_pool{n_threads};
                std::vector<std::future<void>> futures;

                size_t size = 0;
                for (auto& q : queries)
                {
                    futures.emplace_back(pool.execute(this->search, q));
                    size++; // count as generic ranges don't support .size()
                }

                std::vector<result_t> output;
                output.reserve(size);

                for (auto& f : futures)
                    output.push_back(f.get());

                return output;
            }
    };

    // convenient creation function that only takes the ks and picks everything else on it's own
    template<size_t... ks, std::ranges::range text_t>
    auto make_kmer_index(text_t && text, size_t n_threads = 1)
    {
        assert(n_threads > 0);

        using alphabet_t = seqan3::innermost_value_type_t<text_t>;
        using position_t = uint32_t;

        return kmer_index<alphabet_t, position_t, ks...>(std::forward<text_t>(text), n_threads);
    }

}// end of namespace kmer