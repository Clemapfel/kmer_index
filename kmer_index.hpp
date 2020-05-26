//
// Created by Clemens Cords on 2/6/20.
// Copyright (c) 2020 Clemens Cords. All rights reserved.
//

#pragma once

#include <seqan3/alphabet/concept.hpp>
#include <seqan3/alphabet/hash.hpp>
#include <seqan3/range/views/kmer_hash.hpp>
#include <seqan3/range/views/slice.hpp>
#include <seqan3/search/configuration/max_error.hpp>
#include <seqan3/core/type_traits/range.hpp>
#include <seqan3/core/debug_stream.hpp>

#include <type_traits>
#include <cstdint>
#include <cassert>
#include <cmath>
#include <algorithm>
#include <unordered_set>
#include <limits>
#include <thread>
#include <mutex>
#include <unordered_map>

#include <robin_hood.h>

#include "kmer_index_result.hpp"
#include "thread_pool.hpp"

namespace detail
{
    // optimized consteval pow
    constexpr size_t fast_pow(size_t base, size_t exp)
    {
        int result = 1ul;
        for (;;)
        {
            if (exp & 1ul)
                result *= base;
            exp >>= 1ul;
            if (!exp)
                break;
            base *= base;
        }

        return result;

        // reference: https://stackoverflow.com/questions/101439/the-most-efficient-way-to-implement-an-integer-based-power-function-powint-int
    }
}

// represents a kmer index for a single k
template<seqan3::alphabet alphabet_t, size_t k, typename position_t>
class kmer_index_element
{
    static_assert(k > 1, "please specify a valid k");

    private:
        // typedefs for readability
        using result_t = detail::kmer_index_result<alphabet_t, k, position_t>;
        constexpr static size_t _sigma = seqan3::alphabet_size<alphabet_t>;

        // data
        robin_hood::unordered_map<size_t, std::vector<position_t>> _data;

        // needed for subk search edge case
        std::vector<alphabet_t> _last_kmer;
        std::vector<std::vector<position_t>> _last_kmer_refs;

        // hash a query of arbitrary length
        /*
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

        // hash a query of length k, optimized by compiler unwrapping fold expression
        template<typename iterator_t>
        size_t hash_aux_aux(iterator_t query_it, size_t i) const
        {
            return seqan3::to_rank(*query_it) * detail::fast_pow(_sigma, k - i - 1);
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

        const std::vector<position_t>* at(size_t hash) const
        {
            auto it = _data.find(hash);

            if (it != _data.end())
                return &it->second;
            else
                return nullptr;
        }

        template<typename iterator_t>
        void check_last_kmer(iterator_t subk_begin, size_t size, std::vector<const std::vector<position_t>*>& to_fill) const
        {
            // check edge case at last kmer
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
                {
                    to_fill.push_back(&_last_kmer_refs.at(i));  // returns correct position at text.size() - k + i
                }
            }
        }

        template<typename iterator_t>
        std::vector<const std::vector<position_t>*> get_position_for_all_kmer_with_prefix(iterator_t prefix_begin, size_t size) const
        {
            auto it = prefix_begin;
            size_t prefix_hash = 0;
            for (size_t i = 0; i < size; ++i)
                prefix_hash += seqan3::to_rank(*it++) * detail::fast_pow(_sigma, k - i -1);

            size_t lower_bound = 0 + prefix_hash;
            size_t upper_bound = detail::fast_pow(_sigma, size) - (1 / _sigma) + prefix_hash; //c.f. addendum
            size_t step_size = 1;

            std::vector<const std::vector<position_t>*> output;

            // bc stepsize is 1, just add number of possible hashes to lower_bound;
            for (size_t hash = lower_bound; hash < lower_bound + detail::fast_pow(_sigma, k - size); hash += step_size)
            {
                const auto* pos = at(hash);
                if (pos != nullptr)
                    output.push_back(pos);
            }

            check_last_kmer(prefix_begin, size, output);
            return output;

            // Addendum:
            // To calculate set H of hashes for all kmers with suffix s
            // i)   hash of a kmer q_0 q_1 q_2 ... q_k = sum {i=0 to k} rank(q_i) * sigma^(k - i - 1)
            //
            // ii)  for a suffix s_0,..., s_m , the first part of the sum is known:
            //      h_p = sum {i= 0 to m} rank(s_i) * sigma^(k - i - 1)
            //
            // iii) for arbitrary suffix p_0, ..., p_(k-m-1) we observe:
            //      min(H) = hs
            //      by setting all p_i so that rank(p_i) = 0
            //
            // iv)  max(H) = hs + (hash of suffix with all p_i so that r(p_i) = maximum = sigma-1)
            //             = hs + sum {i=k-m to i=k} (sigma-1) * sigma^(k-i-1)
            //             = hs + sigma^m - 1/sigma
            //
            // v)   for two hashes from H h_1 and h_2 so that h_1 < h_2 it is true that
            //      The minimum increase is achieved by setting the last possible char of the suffix of h_2
            //      1 rank higher than the equivalent position of h_1. Because the input of that char to the
            //      value of the hash is 1 * sigma^(k - k-1 - 1) = 1 we know that for all pairs off h_1, h_2 so
            //      the difference h_2 - h_1 is minimal, it holds that h_2 - h_1 = 1
            //
            // vi)  thus to generate all hashes we start with min(H) = hs and stepwise add 1
            //      until we reach max(H) = hs + sigma^m - 1/sigma
        }

    public:

        result_t search(std::vector<alphabet_t>& query) const
        {
            if (query.size() == k)
            {
                const auto* pos = at(hash(query.begin()));
                if (pos)
                    return result_t(pos, this, true);
                else
                    return result_t(this, true);
            }
            else if (query.size() > k)
            {
                size_t rest_n = query.size() % k;

                // get positions for nk parts
                std::vector<const std::vector<position_t>*> nk_positions;

                for (size_t i = 0; i < query.size() - rest_n; i += k)
                {
                    const auto* pos = at(hash(query.begin() + i));

                    if (pos)
                        nk_positions.push_back(pos);
                    else
                        return result_t(this);
                }

                // precompute rest positions
                auto usable = detail::compressed_bitset(nk_positions.back()->size(), true);

                //3748
                if (rest_n > 0)
                {
                    auto rest_results = get_position_for_all_kmer_with_prefix(query.end() - rest_n, rest_n);

                    //auto rest_results_vec = result_t(rest_results, this, true).to_vector();

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
                        {
                            usable.set_1(i);
                            //seqan3::debug_stream << nk_positions.back()->at(i) << " : " << 1 << "\n";
                        }
                        else
                        {
                            usable.set_0(i);
                            //seqan3::debug_stream << nk_positions.back()->at(i) << " : " << 0 << "\n";
                        }

                        ++i;
                    }
                }

                // query.size % k != 0 and query.size < 2*k
                if (nk_positions.size() == 1)
                {
                    result_t output(nk_positions.at(0), this, false);
                    for (size_t i = 0; i < nk_positions.at(0)->size(); ++i)
                        if (usable.at(i))
                            output.should_use(i);

                    return output;
                }

                if (rest_n == 0)
                {
                    result_t output(nk_positions.at(0), this, true);

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
                else
                {
                    result_t output(nk_positions.at(0), this, false);

                    for (size_t start_pos_i = 0; start_pos_i < nk_positions.at(0)->size(); ++start_pos_i)
                    {
                        size_t previous_pos = nk_positions.front()->at(start_pos_i);

                        if (nk_positions.size() > 1)
                        {
                            bool interrupted = false;
                            for (size_t next_pos_i = 1; next_pos_i < nk_positions.size(); ++next_pos_i)
                            {
                                const auto* current = nk_positions.back();
                                auto it = std::find(current->begin(), current->end(), previous_pos += k);

                                if (it == current->end())
                                {
                                    interrupted = true;
                                    break;
                                }

                                // for last k part, also check rest pos
                                if (next_pos_i == nk_positions.size()-1)
                                {
                                    if (not usable.at(it - current->begin()))
                                    {
                                        interrupted = true;
                                        break;
                                    }
                                    else
                                    {
                                        auto vec = usable.to_vector();
                                        break;
                                    }
                                }
                            }

                            if (interrupted)
                                output.should_not_use(start_pos_i);
                            else
                                output.should_use(start_pos_i);
                        }
                    }

                    return output;
                }
            }

            else //query.size() < k
            {
                return result_t(get_position_for_all_kmer_with_prefix(query.begin(), query.size()), this, true);
            }
        }

        // ctors
        ~kmer_index_element() = default;

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

        template<std::ranges::range text_t>
        void create(text_t & text)
        {
            auto hashes = text | seqan3::views::kmer_hash(seqan3::shape{seqan3::ungapped{k}});

            position_t i = 0;
            for (auto h : hashes)
            {
                if (_data.find(h) == _data.end())
                    _data.emplace(h, std::vector<position_t>());

                _data[h].push_back(i);
                ++i;
            }

            size_t text_size = i + k - 1;

            // modify last_kmer_refs with now available text size
            for (auto& ref : _last_kmer_refs)
                ref[0] += text_size - k;

        }

        /*
        // returns bool so fold expression in kmer_index can shortcircuit (c.f. kmer_index exact search comments)
        bool search_if(bool expression, std::vector<alphabet_t> query, std::vector<position_t> &result) const
        {
            if (expression == true)
            {
                result = search(query);
                return true;
            }
            else
                return false;
        }
         */
};

/*

// "multi" kmer index, specifying more than one k depending on output can drastically increase performance [1]
template<seqan3::alphabet alphabet_t, typename position_t, size_t... ks>
class kmer_index
        : protected detail::kmer_index_element<alphabet_t, ks, position_t>...
{
    private:
        // shorthand typedef for readability
        template<size_t k>
        using index_element = detail::kmer_index_element<alphabet_t, k, position_t>;

        // all k used in template, used to determine which element to invoke during query searching
        inline static auto _all_ks = std::vector<size_t>{ks...};
        inline static auto _max_k = *(std::max_element(_all_ks.begin(), _all_ks.end()));

        // search query <= k with index element with optimal k
        std::vector<position_t> search_query_length_subk_or_k(std::vector<alphabet_t> query) const
        {
            assert(query.size() > 0);

            size_t optimal_k = 0;
            for (auto k : _all_ks)
                if (std::labs(k - query.size()) < std::labs(k - optimal_k)) // [2]
                    optimal_k = k;

            std::vector<position_t> result;

            // use boolean fold expression to short-circuit execution and save a few search calls
            (... || index_element<ks>::search_if(ks == optimal_k, query, result));

            return result;
        }

    public:
        // ctor
        template<std::ranges::range text_t>
        kmer_index(text_t && text, size_t n_threads = std::thread::hardware_concurrency(), bool use_hashtable = true)
            : detail::kmer_index_element<alphabet_t, ks, position_t>(use_hashtable)...
        {
            // catch std::hardware_concurrency failing to compute
            if (n_threads == 0)
                n_threads = 1;

            auto pool = thread_pool{n_threads};

            std::vector<std::future<void>> futures;

            // use multiple threads to build index elements at the same time
            (futures.emplace_back(pool.execute(&index_element<ks>::template create<text_t>, static_cast<index_element<ks> *>(this),
                                               std::ref(text))), ...);

            // wait to finish
            for (auto &f : futures)
                f.get();
        }

        // caluclate size (for debugging)
        unsigned long long calculate_size() const
        {
            unsigned long long size = 0;

            size += (index_element<ks>::calculate_size() + ...);
            size += sizeof(_all_ks);

            return size;
        }

        // exact search
        std::vector<position_t> search(std::vector<alphabet_t> query) const
        {
            // if there's only one element, skip choosing which to use
            if (_all_ks.size() == 1)
                return (index_element<ks>::search(query), ...); // expands to single call

            //  if search can be done with exactly one search call
            if (query.size() <= _max_k)
                return search_query_length_subk_or_k(query);

            // else, split query and search each part with corresponding optimal index element [1]
            else
            {
                // if one kmer_index has a way to call search_nk
                std::vector<position_t> result;
                bool possible = (... || index_element<ks>::search_if(query.size() % ks == 0, query, result));

                if (possible)
                    return result;

                // find k so that has optimal conditions to search [1]
                size_t optimal_k = query.size();

                for (auto k : _all_ks)
                    if (std::labs(query.size() % k - k) < std::labs(query.size() % optimal_k - optimal_k))
                        optimal_k = k;

                auto rest = query.size() % optimal_k;
                std::vector<alphabet_t> part_nk{query.begin(), query.begin() + query.size() - rest};
                std::vector<alphabet_t> part_subk{query.begin() + query.size() - rest, query.end()};

                std::vector<position_t> nk_results{};
                (... || index_element<ks>::search_if(ks == optimal_k, part_nk, nk_results));

                std::vector<position_t> subk_results = search_query_length_subk_or_k(part_subk);

                // merge results for first n*k parts and <k rest
                std::vector<position_t> confirmed_positions{};
                for (auto& pos : nk_results)
                {
                    if (std::find(subk_results.begin(), subk_results.end(), pos + query.size() - rest) != subk_results.end())
                        confirmed_positions.push_back(pos);
                }

                return confirmed_positions;

                // [1] search runtime is inversely proportional to r = abs(query.size() - k) c.f. index_element::search
            }
        }

        // serach multiple queries in paralell
        std::vector<std::vector<position_t>> search(
                std::vector<std::vector<alphabet_t>> queries,
                size_t n_threads = std::thread::hardware_concurrency()) const
        {
            using namespace seqan3;

            if (n_threads == 0)
                n_threads = 1;

            auto pool = thread_pool{n_threads};
            std::vector<std::future<std::vector<position_t>>> futures;

            for (auto& q : queries)
                futures.emplace_back(pool.execute(this->search, q));

            std::vector<std::vector<position_t>> output;
            for (auto& f : futures)
                output.emplace_back(f.get());

            return output;
        }
};

// convenient make function that picks template params for you
template<size_t... ks, std::ranges::range text_t>
auto make_kmer_index(text_t&& text)
{
    assert(text.size() < UINT32_MAX && "your text is too large for this configuration, please specify template parameter position_t = uint64_t manually");

    using alphabet_t = seqan3::innermost_value_type_t<text_t>;
    using position_t = uint32_t;

    // generic size since not all ranges support .size()
    size_t size = 0;
    for (const auto _ : text)
        size++;

    // determine if hashtable should be used over map
    bool use_hashtable = false; //text < 10000 or not k > 10; //TODO:: use estimate function ;

    return kmer_index<alphabet_t, position_t, ks...>{
            std::forward<text_t>(text),
            std::thread::hardware_concurrency(),
            use_hashtable};
}

namespace debug
{
    template<size_t... ks, std::ranges::range text_t>
    auto make_kmer_index(text_t&& text, bool use_hashtable, size_t n_threads)
    {
        assert(text.size() < UINT32_MAX && "your text is too large for this configuration, please specify template parameter position_t = uint64_t manually");

        using alphabet_t = seqan3::innermost_value_type_t<text_t>;
        using position_t = uint32_t;

        // generic size since not all ranges support .size()
        size_t size = 0;
        for (const auto _ : text)
            size++;

        return kmer_index<alphabet_t, position_t, ks...>{std::forward<text_t>(text), n_threads, use_hashtable};
    }
}

template<auto first, auto... ks, std::ranges::range text_t>
auto make_kmer_index(text_t text, size_t n_threads = 1)
{
assert(text.size() < UINT32_MAX && "your text is too large for this configuration, please specify template parameter position_t = uint64_t manually");
assert(n_threads > 0);

using alphabet_t = seqan3::innermost_value_type_t<text_t>;
using position_t = uint32_t;

return std::conditional_t<
        std::is_same_v<decltype(first), bool>,
        kmer_index<alphabet_t, position_t, first, ks...>,
        kmer_index<alphabet_t, position_t, true, first, ks...>
>{text, n_threads};
}
 */




