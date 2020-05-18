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

#include <robin_hood.h>

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
#include "thread_pool.hpp"

namespace detail {
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

    public:
        robin_hood::unordered_map<size_t, std::vector<position_t>> _data;

        std::vector<alphabet_t> _first_kmer; // needed for subk search edge case

        constexpr static size_t _sigma = seqan3::alphabet_size<alphabet_t>;

        /*
        template<typename iterator_t>
        static size_t hash(iterator_t query_it)
        {
            size_t hash = 0;
            for (size_t i = 0; i < k; ++i)
                hash += (seqan3::to_rank(*query_it++) * detail::fast_pow(_sigma, k - i - 1));

            return hash;
        }*/

        // hash a query, compiler partially unwraps for loop at compile time
        template<typename iterator_t>
        static size_t hash_aux_aux(iterator_t query_it, size_t i)
        {
            return seqan3::to_rank(*query_it) * detail::fast_pow(_sigma, k - i - 1);
        }

        template<typename iterator_t, size_t... is>
        static size_t hash_aux(iterator_t query_it, std::index_sequence<is...> sequence)
        {
            return (... + hash_aux_aux(query_it++, is));
        }

        template<typename iterator_t>
        static size_t hash(iterator_t query_it)
        {
            return hash_aux(query_it, std::make_index_sequence<k>());
        }

    //protected:

        const std::vector<position_t> _empty;
        std::vector<std::vector<position_t>> _first_kmer_refs;

        const std::vector<position_t>* at(size_t hash) const
        {
            auto it = _data.find(hash);

            if (it != _data.end())
                return &it->second;
            else
                return &_empty;
        }


        template<typename iterator_t>
        std::vector<const std::vector<position_t>*> search_subk(iterator_t suffix_begin, size_t size) const //TODO: iterator
        {
            // generate hashes for all kmers that have query as suffix
            std::vector<std::array<int8_t, k>> rank_sums;
            rank_sums.reserve(detail::fast_pow(_sigma, k - size));

            std::array<int8_t, k> primer;
            auto it = suffix_begin;
            for (size_t i = 0; i < k; ++i)
            {
                if (i >= k - size)
                {
                    primer[i] = seqan3::to_rank(*it++);
                }
                else
                    primer[i] = -1;
            }

            rank_sums.push_back(primer);

            std::array<int8_t, _sigma> current_n_appended{};

            for (int8_t i = k - size - 1; i >= 0; --i)
            {
                // duplicate for next level
                for (auto& sum : rank_sums)
                {
                    for (size_t j = 0; j < _sigma - 1; ++j)
                    {
                        rank_sums.push_back(sum);
                    }
                }

                current_n_appended.fill(0);
                // append in front
                for (auto& sum : rank_sums)
                {
                    sum[i] = current_n_appended[sum.at(i + 1)];
                    current_n_appended[sum.at(i + 1)] += 1;
                }
            }

            std::vector<const std::vector<position_t>*> output;
            output.reserve(rank_sums.size());

            for (auto& h : rank_sums)
            {
                size_t hash = 0;
                for (size_t i = 0; i < k; ++i)
                    hash += (h.at(i) * detail::fast_pow(_sigma, k - i - 1));

                output.push_back(at(hash));
            }

            // check edge case at first kmer
            std::vector<size_t> pos;
            for (size_t i = 0; i < k - size; ++i)
            {
                bool equal = true;
                it = suffix_begin;

                for (size_t j = i; j < i + size; ++j)
                {
                    if (_first_kmer.at(j) != *it)
                    {
                        equal = false;
                        break;
                    }
                    it++;
                }

                if (equal)
                {
                    output.push_back(&_first_kmer_refs.at(i));
                }
            }

            return output;
        }

    public:
        // search any query
        std::vector<position_t> search(std::vector<alphabet_t> query) const
        {
            if (query.size() == k)
            {
                //return std::vector<const std::vector<position_t>*>{at(hash(query))};

                std::vector<position_t> output;
                for (const auto pos : *at(hash(query.begin())))
                    output.push_back(pos);
                return output;
            }

            else if (query.size() > k)
            {
                size_t rest_n = query.size() % k;

                std::vector<const std::vector<position_t>*> positions;

                // positions for first n parts of length k
                for (size_t i = 0; i < query.size() - rest_n; i += k)
                {
                    auto my_hash = hash(query.begin() + i);

                    auto it = _data.find(hash(query.begin() + i));
                    if (it == _data.end())
                        return std::vector<position_t>();
                    else
                        positions.push_back(&(it->second));
                }

                std::vector<const std::vector<position_t>*> rest_positions;

                // position for last part <k
                if (rest_n > 0)
                   rest_positions = search_subk(query.end() - rest_n, rest_n);

                // find out if pos for sections match
                std::vector<position_t> confirmed_positions{};
                for (position_t start_pos : *positions.at(0))
                {
                    position_t previous_pos = start_pos;

                    for (size_t i = 1; i <= positions.size(); ++i)
                    {
                        if (i == positions.size())
                        {
                            if (rest_n > 0)
                            {
                                for (const auto* vec : rest_positions)
                                    if (std::find(vec->begin(), vec->end(), previous_pos + k) !=
                                        vec->end())
                                            confirmed_positions.push_back(start_pos);
                            }
                            else
                                confirmed_positions.push_back(start_pos);

                            break;
                        }

                        const auto* current = positions.at(i);

                        if (std::find(current->begin(), current->end(), previous_pos + k) != current->end())
                            previous_pos += k;
                        else
                            break;
                    }
                }

                return confirmed_positions;
            }

            else //query.size() < k
            {
                //return search_subk();

                std::vector<position_t> output;
                for (const auto* vec : search_subk(query.begin(), query.size()))
                    for (const auto pos : *vec)
                        output.push_back(pos);
                return output;
            }
        }

        // ctors
        ~kmer_index_element() = default;

        template<std::ranges::range text_t>
        kmer_index_element(text_t& text)
                :  _first_kmer(text.begin(), text.begin() + k)
        {
            for (position_t i = 0; i < _first_kmer.size(); ++i)
                _first_kmer_refs.emplace_back(std::vector<position_t>{i});

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
                    _data.emplace(h, std::vector<position_t>({i}));
                else
                    _data[h].push_back(i);

                ++i;
            }
        }

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




