//
// Created by Clemens Cords on 2/6/20.
// Copyright (c) 2020 Clemens Cords. All rights reserved.
//

#pragma once

#include <seqan3/alphabet/concept.hpp>
#include <seqan3/alphabet/hash.hpp>
#include <seqan3/range/views/kmer_hash.hpp>
#include <seqan3/core/type_traits/range.hpp>

#include <type_traits>
#include <cstdint>
#include <cassert>

namespace detail
{
// constexpr pow for longs
constexpr unsigned long long pow_ul(unsigned int base, unsigned int exponent)
{
    unsigned long result = 1;
    while (exponent > 0)
    {
        if (exponent & 1)   // odd exponent
            result = result * base;

        exponent = exponent >> 1;
        base = base * base;
    }

    return result;
}

// picks the smallest datatype that still fits n_bits
template<uint64_t n_bits>
using minimal_memory_t = std::conditional_t<(n_bits < UINT8_MAX), uint8_t,
        std::conditional_t<(n_bits < UINT16_MAX), uint16_t,
                std::conditional_t<(n_bits < UINT32_MAX), uint32_t,
                        uint64_t
                >
        >
>;

// represents a kmer index for a single k
template<seqan3::alphabet alphabet_t,
        size_t k,
        typename hash_t,
        typename position_t,
        bool use_direct_addressing>
class kmer_index_element
{
    static_assert(k > 0 and detail::pow_ul(seqan3::alphabet_size<alphabet_t>, k) < UINT64_MAX,
        "k and alphabet combination are invalid, please choose a smaller non-zero k");

    private:
        // regular mode: minimal memory used, constant runtime for search
        std::unordered_map<hash_t, std::vector<position_t>> _data;

        // direct addressing mode: maximum memory used but o(1) runtime
        std::array<std::vector<position_t>,
                static_cast<int>(use_direct_addressing) * detail::pow_ul(seqan3::alphabet_size<alphabet_t>, k)> _da_data;

        std::vector<alphabet_t> _first_kmer; // needed for subk search edge case

    protected:
        static hash_t hash(std::vector<alphabet_t> query)
        {
            assert(query.size() == k);

            auto shape = seqan3::shape{seqan3::ungapped{static_cast<uint8_t>(query.size())}};
            return static_cast<hash_t>( *(query | seqan3::views::kmer_hash(shape)).begin());
        }

        // split query into parts of length k with rest at the end
        static std::vector<std::vector<alphabet_t>> split_query(std::vector<alphabet_t> query)
        {
            std::vector<std::vector<alphabet_t>> parts{};
            parts.reserve(query.size() / k + 1);

            if (query.size() < k)
                return std::vector{query};

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

            // average lookup o(_data.size() / 2)
            if (not use_direct_addressing)
            {
                auto it = _data.find(hash(query));
                if (it != _data.end())
                    return it->second;
                else
                    return std::vector<position_t>{};
            }
                // lookup always o(1)
            else
                return _da_data[hash(query)];
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
            bool not_equal = false;
            for (size_t i = 0; i < query.size(); ++i)
            {
                if (query.at(i) != _first_kmer.at(i))
                {
                    not_equal = true;
                    break;
                }
            }

            if (not not_equal)
                confirmed_positions.push_back(0);

            auto all_kmer = get_all_kmer_with_suffix(query);
            position_t offset = k - query.size();

            for (auto kmer : all_kmer)
            {
                for (auto &pos : search_query_length_k(kmer))
                    confirmed_positions.push_back(pos + offset);
            }

            return confirmed_positions;
        }

        template<std::ranges::range text_t>
        void construct(text_t &&text)
        {
            auto hashes = text | seqan3::views::kmer_hash(seqan3::shape{seqan3::ungapped{k}});
            _first_kmer = std::vector<alphabet_t>{text.begin() + (text.size() - k), text.end()};

            if (use_direct_addressing)
                _da_data.fill(std::vector<position_t>{});

            position_t pos = 0;
            for (const hash_t hash : hashes)
            {
                if (use_direct_addressing)
                    _da_data[hash].push_back(pos);
                else
                    _data[hash].push_back(pos);

                pos += 1;
            }
        }

        // ctors
        kmer_index_element() = delete;
        ~kmer_index_element() = default;

        template<std::ranges::range text_t>
        kmer_index_element(text_t &&text)
        {
            construct(std::forward<text_t>(text));
        }

        // search any query
        std::vector<position_t> search(std::vector<alphabet_t> query) const
        {
            std::vector<position_t> result;

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

template<seqan3::alphabet alphabet_t,
        typename hash_t,
        typename position_t,
        bool should_use_direct_addressing,
        size_t... ks>
class kmer_index
        : detail::kmer_index_element<alphabet_t, ks, hash_t, position_t, should_use_direct_addressing>...
{
    private:
        template<size_t k>
        using index = detail::kmer_index_element<alphabet_t, k, hash_t, position_t, should_use_direct_addressing>;
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

    public:
        // ctor
        template<std::ranges::range text_t>
        kmer_index(text_t && text)
            : detail::kmer_index_element<alphabet_t, ks, hash_t, position_t, should_use_direct_addressing>(text)...
        {
        }

        // exact search
        std::vector<position_t> search(std::vector<alphabet_t> query) const
        {
            if (_all_ks.size() == 1)
                return (index<ks>::search(query), ...); // expands to one call

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

                // if not, split optimally so that as many segments as possible can be searched as kmers
                size_t optimal_k = query.size();
                for (auto k : _all_ks)
                    if (optimal_k % query.size() < k % query.size()) //[1]
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

                // [1] : kmer runtime is inversely proprotional to abs(query.size() - k)
                // so k that's closest to query.size() should be chosen
            }
        }

    TODO: approximate search
};

// MAKE FUNCTION for easier use by defaulting all the template params
// can't be defaulted in the class template because param pack has to be last and you can't use default
// template params along with param pack

template<size_t... ks, bool use_direct_addressing, typename hash_t, typename position_t, std::ranges::range text_t>
auto make_kmer_index(text_t text)
{
    return kmer_index<seqan3::innermost_value_type_t<text_t>, hash_t, position_t, use_direct_addressing, ks...>{text};
}

template<size_t... ks, bool use_direct_addressing, std::ranges::range text_t>
auto make_kmer_index(text_t text)
{
    assert(text.size() < UINT32_MAX && "please specify template parameter position_t = uint64_t manually");

    using alphabet_t = seqan3::innermost_value_type_t<text_t>;
    using hash_t = detail::minimal_memory_t<std::max(detail::pow_ul(ks, seqan3::alphabet_size<alphabet_t>)..., 0ull)>;
    using position_t = uint32_t;

    return kmer_index<alphabet_t, hash_t, position_t, use_direct_addressing, ks...>{text};
}

template<size_t... ks, std::ranges::range text_t>
auto make_kmer_index(text_t text)
{
    assert(text.size() < UINT32_MAX && "please specify template parameter position_t = uint64_t manually");

    // code copy pasted to prevent circular type deduction (not possible: make_kmer_index<ks..., false>(text))
    using alphabet_t = seqan3::innermost_value_type_t<text_t>;
    using hash_t = detail::minimal_memory_t<std::max(detail::pow_ul(ks, seqan3::alphabet_size<alphabet_t>)..., 0ull)>;
    using position_t = uint32_t;

    return kmer_index<alphabet_t, hash_t, position_t, false, ks...>{text};
}



