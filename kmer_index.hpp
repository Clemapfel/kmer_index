//
// Created by Clemens Cords on 2/6/20.
// Copyright (c) 2020 Clemens Cords. All rights reserved.
//

#pragma once

#include <seqan3/alphabet/concept.hpp>
#include <seqan3/alphabet/hash.hpp>
#include <seqan3/range/views/kmer_hash.hpp>
#include <seqan3/core/type_traits/range.hpp>
#include <seqan3/core/debug_stream.hpp>

#include <cmath>
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
        // represents a kmer-index for a single set k
        // alphabet_t   :   the alphabet of the text
        // position_t   :   the primitive used for positional indices
        // k            :   the k user for hashing kmers
        template<seqan3::alphabet alphabet_t, typename position_t, size_t k>
        class kmer_index_element
        {
            static_assert(k > 0 and k < 64 / log2(seqan3::alphabet_size<alphabet_t>),
                    "the hashspace for the current k cannot be represented with only a 64-bit integer. Please specify a valid k");

            friend class kmer_index;

            private:
                // typedefs for readabilty
                using result_t = kmer_index_result<position_t>;
                constexpr static size_t _sigma = seqan3::alphabet_size<alphabet_t>;

                robin_hood::unordered_map<size_t, std::vector<position_t>> _data;

                // hash a query of length k
                // optimization through contesxpr unwrapping parts of fold expression
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

                // access data
                const std::vector<position_t>* at(size_t hash) const
                {
                    auto it = _data.find(hash);

                    if (it != _data.end())
                        return &it->second;
                    else
                        return nullptr;
                }

                // check last kmer to account for edge case
                std::vector<alphabet_t> _last_kmer;
                std::vector<std::vector<position_t>> _last_kmer_refs;

                template<typename iterator_t>
                void check_last_kmer(iterator_t subk_begin, size_t size,
                                     std::vector<const std::vector<position_t>*>& to_fill) const
                {
                    for (size_t i = 1; i < k - size + 1; ++i)
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
                    }
                }

                // access positions based on prefix of length < k
                template<typename iterator_t>
                std::vector<const std::vector<position_t>*>
                get_position_for_all_kmer_with_prefix(iterator_t prefix_begin, size_t size) const
                {
                    if (fast_pow(_sigma, k-size) > 1e7)
                    {
                        throw std::invalid_argument("query size too low for specified k");
                    }

                    auto it = prefix_begin;
                    size_t prefix_hash = 0;
                    for (size_t i = 0; i < size; ++i)
                    {
                        prefix_hash += seqan3::to_rank(*it++) * fast_pow(_sigma, k - i - 1);
                    }

                    size_t lower_bound = 0 + prefix_hash;
                    size_t n_hashes =  kmer::detail::fast_pow(_sigma, k - size);
                    size_t upper_bound = lower_bound + n_hashes;
                    size_t step_size = 1;

                    std::vector<const std::vector<position_t>*> output;

                    for (size_t hash = lower_bound;
                         hash < upper_bound; hash += step_size)
                    {
                        const auto* pos = at(hash);
                        if (pos != nullptr)
                            output.push_back(pos);
                    }

                    check_last_kmer(prefix_begin, size, output);
                    return output;
                }

            protected:
                // CTOR protected because user should only engage with kmer_index_element<k> through kmer_index<k>
                kmer_index_element() = default;

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

                    assert(i + k -1 < std::numeric_limits<position_t>::max() &&
                        "your text is too large for this configuration");

                    position_t text_size = position_t(i) + k - 1;

                    _last_kmer = std::vector<alphabet_t>(text.end() - k, text.end());
                    _last_kmer_refs = std::vector<std::vector<position_t>>();

                    for (position_t j = 0; j < _last_kmer.size(); ++j)
                        _last_kmer_refs.emplace_back(std::vector<position_t>{j + text_size - position_t(k)});
                }

            public:
                template<typename iterator_t>
                const std::vector<position_t>* search_k(iterator_t it) const
                {
                    const auto *pos = at(hash(it));
                    if (pos)
                        return pos;
                    else
                        return nullptr;
                }

                // search any query
                virtual result_t search(std::vector<alphabet_t>& query) const
                {
                    assert(query.size() > 0);

                    // query size exactly k
                    if (query.size() == k)
                    {
                        const auto* pos = at(hash(query.begin()));
                        if (pos)
                            return result_t(pos, true, BYPASS_BITMASK::YES);
                        else
                            return result_t();
                    }
                    // query size m > k
                    else if (query.size() > k)
                    {
                        size_t rest_n = query.size() % k;

                        // get positions for nk parts
                        std::vector<const std::vector<position_t>*> nk_positions;

                        int last_hash = -1;

                        for (size_t i = 0; i < query.size() - rest_n; i += k)
                        {
                            size_t h = hash(query.begin() + i);
                            const auto* pos = (h != last_hash ? at(h) : nk_positions.back());

                            if (pos)
                                nk_positions.push_back(pos);
                            else
                                return result_t();

                            last_hash = h;
                        }

                        // get positions for rest
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
                            result_t output(nk_positions.at(0), false, BYPASS_BITMASK::NO);
                            for (size_t i = 0; i < nk_positions.at(0)->size(); ++i)
                                if (usable.at(i))
                                    output.should_use(i);

                            return output;
                        }

                        // query.size % k == 0
                        else if (rest_n == 0)
                        {
                            result_t output(nk_positions.at(0), true, BYPASS_BITMASK::NO);

                            for (size_t start_pos_i = 0; start_pos_i < nk_positions.front()->size(); ++start_pos_i)
                            {
                                size_t previous_pos = nk_positions.front()->at(start_pos_i);
                                bool should_use = true;

                                for (size_t next_pos_i = 1; next_pos_i < nk_positions.size(); ++next_pos_i)
                                {
                                    auto* current = nk_positions.at(next_pos_i);

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
                            result_t output(nk_positions.at(0), true, BYPASS_BITMASK::NO);

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

                    // query.size() < k
                    else
                    {
                        return result_t(get_position_for_all_kmer_with_prefix(query.begin(), query.size()));
                    }
                }
        };
    } // end of namespace detail

    template<seqan3::alphabet alphabet_t, typename position_t, size_t... ks>
    class kmer_index
        : public detail::kmer_index_element<alphabet_t, position_t, ks>...
    {
        private:
            template<size_t k>
            using index_element_t = detail::kmer_index_element<alphabet_t, position_t, k>;
            using index_t = kmer_index<alphabet_t, position_t, ks...>;
            using result_t = detail::kmer_index_result<position_t>;

            inline static std::vector<size_t> _all_ks = std::vector<size_t>{ks...};

            // ###
            // setup such that kmer_index can call a specific search function of one of it's elements
            // by specifying a non-constexpr k during runtime

            // we first fill an array with functions for an individual kmer_index_elements' search (and search_k)
            template<size_t k>
            result_t call_search(std::vector<alphabet_t>& query) const
            {
                return static_cast<const index_element_t<k>*>(this)->index_element_t<k>::search(query);
            }

            using search_fn = result_t(kmer_index<alphabet_t, position_t, ks...>::*)(
                    std::vector<alphabet_t>&) const;

            const std::array<search_fn, sizeof...(ks)> _search_fns = {
                    (&kmer_index<alphabet_t, position_t, ks...>::call_search<ks>)...};

                        template<size_t k>
            const std::vector<position_t>* call_search_k(typename std::vector<alphabet_t>::iterator query_begin) const
            {
                return static_cast<const index_element_t<k>*>(this)->index_element_t<k>::search_k(query_begin);
            }

            using search_k_fn = const std::vector<position_t>*(kmer_index<alphabet_t, position_t, ks...>::*)(
                    typename std::vector<alphabet_t>::iterator) const;

            const std::array<search_k_fn, sizeof...(ks)> _search_k_fns = {
                    (&kmer_index<alphabet_t, position_t, ks...>::call_search_k<ks>)...};

            // then fill an array such that _k_to_search_fns_i[k] returns the correct
            // index to get a function from _search_fns
            // now we can call this->*(_search_fns[k])(query); where k ist non-constexpr and
            // invoke the corresponding kmer_index_element with RVO
            std::array<int, std::max({ks...}) + 1> _k_to_search_fns_i;

            void setup_k_to_search_fn()
            {
                std::fill(_k_to_search_fns_i.begin(), _k_to_search_fns_i.end(), -1);
                size_t i = 0;
                for (size_t k : std::vector<size_t>{ks...})
                {
                    _k_to_search_fns_i[k] = i;
                    ++i;
                }
            }

            // calculate which query lengths to search with which k
            inline static size_t _query_size_range = 10000;
            inline static size_t _max_possible_k = 32;

            std::array<std::vector<size_t>, _query_size_range> _optimal_nk_sum;
            std::array<bool, _query_size_range> _use_multi_search_scheme;

            void choose_search_scheme()
            {
                // TODO: do this at compile time
                std::sort(_all_ks.begin(), _all_ks.end(), [](size_t a, size_t b) -> bool { return a > b; });

                std::vector<size_t> high_ks;
                for (size_t i : _all_ks)
                    if (i >= 9)
                        high_ks.push_back(i);

                std::fill(_optimal_nk_sum.begin(), _optimal_nk_sum.end(), std::vector<size_t>());

                std::fill(_use_multi_search_scheme.begin(), _use_multi_search_scheme.end(), false);

                for (size_t k : high_ks)
                {
                    _optimal_nk_sum[k] = {k};
                    _use_multi_search_scheme[k] = true;
                }

                for (size_t q =_all_ks.front()+1; q < _query_size_range; ++q)
                {
                    bool found = false;
                    for (size_t k : high_ks)
                    {
                        if (not _optimal_nk_sum[q - k].empty())
                        {
                            _optimal_nk_sum[q] = _optimal_nk_sum[q - k];
                            _optimal_nk_sum[q].push_back(k);
                            found = true;
                            break;
                        }
                    }

                    if (found)
                        _use_multi_search_scheme[q] = true;
                }

                for (size_t q = 0; q < _query_size_range; ++q)
                {
                    if (not _optimal_nk_sum[q].empty())
                        continue;

                    if (q < _all_ks.front())
                    {
                        size_t optimal_k = _all_ks.front();
                        for (size_t k : _all_ks)
                        {
                            if (q <= k and (k - q < optimal_k - q))
                            {
                                optimal_k = k;
                                continue;
                            }
                        }
                        _optimal_nk_sum[q] = {optimal_k};
                    }
                    else
                    {
                        size_t optimal_k = _all_ks.front();
                        for (size_t k : _all_ks)
                        {
                            if ((ceil(q / float(k)) * k - q) <
                                (ceil(q / float(optimal_k)) * optimal_k - q))
                                optimal_k = k;
                        }

                        _optimal_nk_sum[q] = {optimal_k};
                    }
                }
            }

        public:
            // CTOR
            template<std::ranges::range text_t>
            kmer_index(text_t& text, size_t n_threads = std::max(std::thread::hardware_concurrency(), 1u))
                    : index_element_t<ks>()...
            {
                auto pool = detail::thread_pool{n_threads};
                std::vector<std::future<void>> futures;
                (futures.emplace_back(
                        pool.execute(&index_element_t<ks>::template create<text_t>,
                                     static_cast<index_element_t<ks>*>(this),
                                     std::ref(text))), ...);
                for (auto& f : futures)
                    f.get();

                setup_k_to_search_fn();
                choose_search_scheme();
            }

            void extend_query_size_range(size_t new_maximum)
            {
                _query_size_range = new_maximum;
                choose_search_scheme();
            }

            // search any query
            result_t search(std::vector<alphabet_t>& query) const
            {
                if (query.size() > _query_size_range)
                    throw(std::invalid_argument("query size exceed the maximum size "
                        + std::to_string(_query_size_range) + " specified"));

                // if rest present just use regular searching
                if (not _use_multi_search_scheme[query.size()] or _all_ks.size() == 1)
                    return (this->*(_search_fns[_k_to_search_fns_i.at(_optimal_nk_sum.at(query.size()).at(0))]))(query);

                // split query into kmers with different k and search each part with appropriate index element
                std::vector<const std::vector<position_t>*> nk_positions;
                size_t last_k = 0;
                for (size_t current_k : _optimal_nk_sum[query.size()])
                {
                    const auto* pos = (this->*_search_k_fns[_k_to_search_fns_i.at(current_k)])(query.begin() + last_k);
                    if (pos)
                        nk_positions.push_back(pos);
                    else
                        return result_t();

                    last_k = current_k;
                }

                if (nk_positions.size() == 1)
                    return result_t(nk_positions.at(0), true, detail::BYPASS_BITMASK::YES);

                result_t output(nk_positions.at(0), true, detail::BYPASS_BITMASK::NO);

                // crossreference the part positions
                size_t nk_sum_i = 0;
                for (size_t start_pos_i = 0; start_pos_i < nk_positions.front()->size(); ++start_pos_i)
                {
                    size_t previous_pos = nk_positions.front()->at(start_pos_i);

                    bool interrupted = false;
                    for (size_t next_pos_i = 1; next_pos_i < nk_positions.size(); ++next_pos_i)
                    {
                        const auto* current = nk_positions.at(next_pos_i);
                        auto it = std::lower_bound(current->begin(), current->end(), previous_pos += _optimal_nk_sum.at(query.size()).at(nk_sum_i));

                        if (*it != previous_pos)
                        {
                            interrupted = true;
                            break;
                        }
                    }

                    if (interrupted)
                        output.should_not_use(start_pos_i);
                }

                return output;
            }

            // overload for rvalue
            result_t search(std::vector<alphabet_t>&& query) const
            {
                auto hold = query;
                search(hold);
            }
    };

    // convenient creation function that only takes the ks and picks everything else on it's own
    template<size_t... ks, std::ranges::range text_t>
    auto make_kmer_index(text_t && text, size_t n_threads = std::thread::hardware_concurrency())
    {
        assert(n_threads > 0);

        using alphabet_t = seqan3::range_innermost_value_t<text_t>;
        using position_t = uint32_t;
        using hash_t = uint64_t;

        return kmer_index<alphabet_t, position_t, ks...>(std::forward<text_t>(text), n_threads);
    }

} // end of namespace kmer