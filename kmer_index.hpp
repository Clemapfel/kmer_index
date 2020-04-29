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
#include "thread_pool.hpp"

namespace detail
{
template<seqan3::alphabet alphabet_t, size_t k>
size_t hash_query(std::vector<alphabet_t> query)
{
    assert(query.size() == k);

    auto hashes = query | seqan3::views::kmer_hash(seqan3::shape{seqan3::ungapped{k}});
    return *(hashes.begin());

    /*

    size_t hash = 0;
    for (size_t i = 0; i < k; ++i)
        hash += seqan3::to_rank(query.at(i)) * pow(k, k-i-1);

    return hash;
     */
}

bool USE_DA_HEURISTIC = true;

// represents a kmer index for a single k
template<seqan3::alphabet alphabet_t, size_t k, typename position_t, bool use_hashtable=true>
class kmer_index_element
{
    private:
        // custom hash table for direct addressing
        class kmer_hash_table
        {
            private:
                std::vector<std::vector<position_t>> _data;

                size_t _text_length;

            protected:
                void print_hashtable() const
                {
                    seqan3::debug_stream << "i" << "\t" << "kmer_hash" << "\t" << "pos" << "\n";

                    size_t i = 0;
                    size_t kmer_hash = _min_hash;
                    for (auto vec : _data)
                    {
                        seqan3::debug_stream << i << "\t" << kmer_hash << "\t" << vec << "\n";
                        i += 1;
                        kmer_hash += 1;
                    }

                    size_t n_empty_slots = 0;
                    for (auto row : _data)
                        if (row.empty())
                            n_empty_slots++;

                    seqan3::debug_stream << "Printed DA Hashtable for k = " << k << ", text_length = " << _text_length << "." << "\n";
                    seqan3::debug_stream << "min hash: " << _min_hash << " , max hash: " << _max_hash
                                         << " (" << ((n_empty_slots / float(_data.size()))) * 100 << "% entries empty)\n";
"\n";
                    seqan3::debug_stream << "______________________________________________\n";
                }

                size_t _min_hash, _max_hash;
                size_t _n_different_hashes;

                // run heuristic to size DA table
                template<std::ranges::range text_t>
                void estimate_hash_range(text_t text, float sample_coverage = 1)
                {
                    if (not USE_DA_HEURISTIC)
                    {
                        for (auto h : text | seqan3::views::kmer_hash(seqan3::shape{seqan3::ungapped{k}}))
                        {
                            _min_hash = h;
                            _max_hash = h;
                            break;
                        }
                        return;
                    }

                    _min_hash = std::numeric_limits<size_t>::max();
                    _max_hash = 0;
                    std::unordered_set<size_t> hashes;

                    for (auto h : text | seqan3::views::kmer_hash(seqan3::shape{seqan3::ungapped{k}}))
                    {
                        if (h < _min_hash)
                            _min_hash = h;
                        if (h > _max_hash)
                            _max_hash = h;

                        hashes.insert(h);
                    }

                    assert(hashes.size() >= 1);
                    _n_different_hashes = hashes.size();
                }

                // construct hash table
                template<std::ranges::range text_t>
                void construct(text_t text)
                {
                    _text_length = text.size();

                    auto hashes = text | seqan3::views::kmer_hash(seqan3::shape{seqan3::ungapped{k}});
                    estimate_hash_range(text);

                    _data.reserve(_max_hash - _min_hash);

                    for (size_t i = 0; i <= (_max_hash - _min_hash); ++i)
                        _data.emplace_back();

                    position_t i = 0;
                    for (auto h : hashes)
                    {
                        // resize if necessary
                        if (h < _min_hash)
                        {
                            _data.insert(_data.begin(), std::labs(h - _min_hash), std::vector<position_t>{});
                            _min_hash = h;
                        }
                        else if (h > _max_hash)
                        {
                            _data.insert(_data.end(), h - _max_hash, std::vector<position_t>{});
                            _max_hash = h;
                        }

                        _data[h - _min_hash].push_back(i);
                        i += 1;
                    }

                    //seqan3::debug_stream << "kmer element for k = " << k << " finished construction.\n";
                }

            public:
                kmer_hash_table() = default;

                // create from text
                template<std::ranges::range text_t>
                kmer_hash_table(text_t text)
                {
                    estimate_hash_range(text);
                    construct(text);
                }

                std::vector<position_t> at(std::vector<alphabet_t> query) const
                {
                    assert(query.size() == k);

                    auto hash = hash_query<alphabet_t, k>(query);

                    if (hash < _min_hash or hash > _max_hash)
                        return std::vector<position_t>{};
                    else
                        return _data.at(hash - _min_hash);
                };
        };

        kmer_hash_table _da_data;
        std::unordered_map<size_t, std::vector<position_t>> _map_data;

        std::vector<alphabet_t> _first_kmer; // needed for subk search edge case

    protected:
        // split query into parts of length k with rest at the end
        static std::vector<std::vector<alphabet_t>> split_query(std::vector<alphabet_t> query)
        {
            std::vector<std::vector<alphabet_t>> parts{};
            parts.reserve(query.size() / k + 1);

            if (query.size() < k)
                return {query};

            for (size_t offset = 0; offset <= query.size() - k; offset += k)
                parts.emplace_back(query.begin() + offset, query.begin() + offset + k);

            if (query.size() % k != 0)
                parts.emplace_back(query.end() - k + 1, query.end());

            return parts;
        }

        // get all possible kmers that have query as suffix, needed for subk search
        static std::unordered_set<std::vector<alphabet_t>> get_all_kmer_with_suffix(std::vector<alphabet_t> sequence)
        {
            assert(sequence.size() <= k);

            std::vector<alphabet_t> all_letters{};

            for (size_t i = 0; i < alphabet_t::alphabet_size; ++i)
                all_letters.push_back(seqan3::assign_rank_to(i, alphabet_t{}));

            std::unordered_set<std::vector<alphabet_t>> output{};

            size_t size = pow(alphabet_t::alphabet_size, (sequence.size()));

            auto current = output;
            current.reserve(size);

            for (auto letter : all_letters)
                output.insert({letter});

            for (size_t i = 1; i < k; i++)
            {
                current.swap(output);
                output.clear();

                if (i < k - sequence.size())
                {
                    for (auto seq : current)
                    {
                        for (auto letter : all_letters)
                        {
                            auto temp = seq;
                            temp.push_back(letter);
                            output.insert(temp);
                        }
                    }
                }
                else
                {
                    for (auto seq : current)
                    {
                        auto temp = seq;
                        temp.push_back(sequence.at(i - (k - sequence.size())));
                        output.insert(temp);
                    }
                }
            }

            return output;
        }

        // exact search for query.size() == k
        std::vector<position_t> search_query_length_k(std::vector<alphabet_t> query) const
        {
            assert(query.size() == k);

            if (use_hashtable)
                return _da_data.at(query);
            else
            {
                auto it = _map_data.find(hash_query<alphabet_t, k>(query));
                if (it != _map_data.end())
                    return it->second;
                else
                    return std::vector<position_t>();
            }
        }

        // exact search for query.size() % k == 0
        std::vector<position_t> search_query_length_nk(std::vector<alphabet_t> query) const
        {
            assert(query.size() % k == 0);

            auto query_parts = split_query(query);

            std::vector<std::vector<position_t>> positions{};
            for (auto &part : query_parts)
            {
                positions.push_back(search_query_length_k(part));

                if (positions.back().empty())
                    return std::vector<position_t>{}; // if one segment is missing, no match
            }

            std::vector<position_t> confirmed_positions{};
            for (auto &start_pos : positions.at(0))
            {
                position_t previous_pos = start_pos;

                for (size_t i = 1; i <= positions.size(); ++i)
                {
                    if (i == positions.size())
                    {
                        confirmed_positions.push_back(start_pos);
                        break;
                    }

                    const auto current = positions.at(i);

                    if (std::find(current.begin(), current.end(), previous_pos + k) != current.end())
                        previous_pos += k;
                    else
                        break;
                }
            }

            return confirmed_positions;
        }

        // exact search for query.size() < k
        std::vector<position_t> search_query_length_subk(std::vector<alphabet_t> query) const
        {
            assert(query.size() < k);

            std::vector<position_t> confirmed_positions{};

            // edge case: query at the very beginning
            for (size_t i = 0; i < query.size(); ++i)
            {
                if (_first_kmer.at(i) == query.at(0))
                {
                    bool equal = true;
                    for (size_t j = 0; j < query.size(); ++j)
                    {
                        if (i + j >= _first_kmer.size() || _first_kmer.at(i + j) != query.at(j))
                        {
                            equal = false;
                            break;
                        }
                    }
                    if (equal)
                        confirmed_positions.push_back(i);
                }
            }

            auto all_kmer = get_all_kmer_with_suffix(query);
            //seqan3::debug_stream << "all kmer with suffix (" << query << ") " << all_kmer << "\n";

            position_t offset = k - query.size();

            for (auto kmer : all_kmer)
            {
                for (auto &pos : search_query_length_k(kmer))
                    confirmed_positions.push_back(pos + offset);
            }

            return confirmed_positions;
        }

        // ctors
        ~kmer_index_element() = default;
        kmer_index_element() = default;

        template<std::ranges::range text_t>
        kmer_index_element(text_t & text)
        {
            create(text);
        }

        template<std::ranges::range text_t>
        void create(text_t & text)
        {
            _first_kmer = std::vector<alphabet_t>(text.begin(), text.begin() + k);

            if (use_hashtable)
            {
                _da_data = kmer_hash_table{text};
            }
            else {
                auto hashes = text | seqan3::views::kmer_hash(seqan3::shape{seqan3::ungapped{k}});

                position_t i = 0;
                for (auto h : hashes) {
                    if (_map_data.find(h) == _map_data.end()) {
                        _map_data.insert(std::make_pair(h, std::vector<position_t>({i})));
                    }
                    else
                        _map_data[h].push_back(i);
                    ++i;
                }

                //seqan3::debug_stream << _map_data << "\n";
            }
        }

        // search any query
        std::vector<position_t> search(std::vector<alphabet_t> query) const
        {
            std::vector<position_t> result;

            if (query.size() % k >= 10)
                seqan3::debug_stream << "[WARNING] searching query (length = " << query.size() << ") with kmer index for k = " << k << " may cause runtime or memory issues. Please choose a k and query combination so that query.size() % k is as high as possible.\n";

            if (query.size() == k)
                return search_query_length_k(query);

            else if (query.size() % k == 0)
                return search_query_length_nk(query);

            else if (query.size() < k)
            {
                // subk results need to be sorted since get_kmers_with_suffix does not return in order
                auto result = search_query_length_subk(query);
                std::sort(result.begin(), result.end());
                return result;
            }
            else
            {
                auto rest = query.size() % k;
                std::vector<alphabet_t> part_nk{query.begin(), query.begin() + query.size() - rest};
                std::vector<alphabet_t> part_subk{query.begin() + query.size() - rest, query.end()};

                auto nk_results = search_query_length_nk(part_nk);
                auto subk_results = search_query_length_subk(part_subk);

                std::vector<position_t> confirmed_positions{};

                for (auto &pos : nk_results)
                {
                    if (std::find(subk_results.begin(), subk_results.end(), pos + query.size() - rest) !=
                        subk_results.end())
                        confirmed_positions.push_back(pos);
                }

                return confirmed_positions;
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
}

template<seqan3::alphabet alphabet_t, typename position_t, bool use_hashtable, size_t... ks>
class kmer_index
        : protected detail::kmer_index_element<alphabet_t, ks, position_t, use_hashtable>...
{
    private:
        template<size_t k>
        using index = detail::kmer_index_element<alphabet_t, k, position_t, use_hashtable>;
        inline static auto _all_ks = std::vector<size_t>{ks...};

        // search query <= k with index element with optimal k
        std::vector<position_t> search_query_length_subk_or_k(std::vector<alphabet_t> query) const
        {
            assert(query.size() > 0);

            size_t optimal_k = 0;
            for (auto k : _all_ks)
                if (std::labs(k - query.size()) < std::labs(k - optimal_k))
                    optimal_k = k;

            std::vector<position_t> result;
            (... || index<ks>::search_if(ks == optimal_k, query, result));
            // use boolean fold expression to short-circuit execution and save a few search calls
            return result;
        }

        inline static size_t _n_threads = 1;

    public:

        // set number of concurrent threads used for search and construction.
        // Defaults to std::thread::hardware_concurrency
        static void set_thread_count(size_t n_threads)
        {
            assert(n_threads > 0);
            _n_threads = n_threads;
        }

        // ctor
        template<std::ranges::range text_t>
        kmer_index(text_t && text, size_t n_threads = std::thread::hardware_concurrency())
            : detail::kmer_index_element<alphabet_t, ks, position_t, use_hashtable>()...
        {

            auto pool = thread_pool{n_threads};

            std::vector<std::future<void>> futures;

            // use multiple threads to build index elements
            (futures.emplace_back(pool.execute(&index<ks>::template create<text_t>, static_cast<index<ks> *>(this),
                                               std::ref(text))), ...);

            // wait to finish
            for (auto &f : futures)
                f.get();
        }

        // exact search
        std::vector<position_t> search(std::vector<alphabet_t> query) const
        {
            if (_all_ks.size() == 1)
                return (index<ks>::search(query), ...); // expands to single call

            size_t max_k = 0;
            for (auto k : _all_ks)
                if (k > max_k)
                    max_k = k;

            if (query.size() <= max_k)
            {
                return search_query_length_subk_or_k(query);
            }
            else
            {
                // if one kmer_index has a way to call search_nk
                std::vector<position_t> result;
                bool possible = (... || index<ks>::search_if(query.size() % ks == 0, query, result));

                if (possible)
                    return result;

                // find k so that rest at the end is as close to that k as possible [1]
                size_t optimal_k = query.size();
                for (auto k : _all_ks)
                    if (std::labs(query.size() % k - k) < std::labs(query.size() % optimal_k - optimal_k))
                        optimal_k = k;

                auto rest = query.size() % optimal_k;
                std::vector<alphabet_t> part_nk{query.begin(), query.begin() + query.size() - rest};
                std::vector<alphabet_t> part_subk{query.begin() + query.size() - rest, query.end()};

                std::vector<position_t> nk_results{};
                (... || index<ks>::search_if(ks == optimal_k, part_nk, nk_results));

                std::vector<position_t> subk_results = search_query_length_subk_or_k(part_subk);

                std::vector<position_t> confirmed_positions{};

                for (auto& pos : nk_results)
                {
                    if (std::find(subk_results.begin(), subk_results.end(), pos + query.size() - rest) != subk_results.end())
                        confirmed_positions.push_back(pos);
                }

                return confirmed_positions;

                // [1] : kmer runtime is inversely proprotional to abs(query.size() - k) bc the bigger the difference,
                // the bigger is the set of kmers that have to be searched as possiblities
                // so k that's closest to query.size() should be chosen
            }
        }

        // serach multiple queries in paralell
        std::vector<std::vector<position_t>> search(std::vector<std::vector<alphabet_t>> queries) const
        {
            using namespace seqan3;

            auto pool = thread_pool{_n_threads};
            std::vector<std::future<std::vector<position_t>>> futures;

            for (auto& q : queries)
                futures.emplace_back(pool.execute(this->search, q));

            // wait to finish
            std::vector<std::vector<position_t>> output;
            for (auto& f : futures)
                output.emplace_back(f.get());

            return output;

            /*

            std::vector<std::vector<position_t>> results;
            results.reserve(queries.size());
            results.assign(queries.size(), std::vector<position_t>());

            std::vector<std::thread> threads;

            size_t i = 0;
            for (auto q : queries)
            {
                auto& current = results[i];
                threads.emplace_back([this, q, &current]()
                    {
                        current = search(q);
                    });
                i++;
            }

            for (auto& thr : threads)
                thr.join();

            return results;
             */
        }

};

template<auto first, auto... ks, std::ranges::range text_t>
auto make_kmer_index(text_t text)
{
    assert(text.size() < UINT32_MAX && "your text is too large for this configuration, please specify template parameter position_t = uint64_t manually");

    using alphabet_t = seqan3::innermost_value_type_t<text_t>;
    using position_t = uint32_t;

  return std::conditional_t<
    std::is_same_v<decltype(first), bool>,
    kmer_index<alphabet_t, position_t, first, ks...>,
    kmer_index<alphabet_t, position_t, true, first, ks...>
  >{text};
}

template<auto first, auto... ks, std::ranges::range text_t>
auto make_kmer_index(text_t text, size_t n_threads)
{
assert(text.size() < UINT32_MAX && "your text is too large for this configuration, please specify template parameter position_t = uint64_t manually");

using alphabet_t = seqan3::innermost_value_type_t<text_t>;
using position_t = uint32_t;

return std::conditional_t<
        std::is_same_v<decltype(first), bool>,
        kmer_index<alphabet_t, position_t, first, ks...>,
        kmer_index<alphabet_t, position_t, true, first, ks...>
>{text, n_threads};
}



