// Copyright (c) 2020 Clemens Cords. All rights reserved.

#pragma once

#include <seqan3/alphabet/concept.hpp>
#include <seqan3/alphabet/hash.hpp>
#include <seqan3/range/views/kmer_hash.hpp>
#include <seqan3/range/views/slice.hpp>
#include <seqan3/search/configuration/max_error.hpp>
#include <seqan3/core/type_traits/range.hpp>
#include <seqan3/core/debug_stream.hpp>

#include <robin_hood.h>

#include <iterator>

namespace detail
{
    // fast pow

}

template<seqan3::alphabet alphabet_t, size_t k>
class minimum_kmer_index
{
    static_assert(k > 1, "please specify a valid k");

    using position_t = uint32_t;

    private:
        robin_hood::unordered_map<size_t, std::vector<position_t>> _data;

        // needed for edge case in sub_k
        std::vector<alphabet_t> _first_kmer;
        const std::vector<position_t> _zero = {0};

        // optimized pow used for hashing
        size_t fast_pow(size_t base, size_t exp)
        {
            int result = 1;
            for (;;)
            {
                if (exp & 1)
                    result *= base;
                exp >>= 1;
                if (!exp)
                    break;
                base *= base;
            }

            return result;

            // reference: https://stackoverflow.com/questions/101439/the-most-efficient-way-to-implement-an-integer-based-power-function-powint-int
        }

        // hash a kmer
        template<typename iterator_t>
        size_t hash(iterator_t& query_it)
        {
            size_t hash = 0;
            for (size_t i = 0; i < k; ++i)
            {
                hash += seqan3::to_rank(*(query_it)) * fast_pow(k, k - i - 1);
                query_it++;
            }

            return hash;
        }


        // find position of all kmer with query as suffix
        template<typename iterator_t>
        std::vector<const std::vector<position_t>*> search_subk(iterator_t query_begin, size_t size)
        {
            iterator_t backup_start_it = query_begin;   // deep copy for later

            // find latter component of hash
            size_t suffix_hash = 0;
            for (size_t i = k - size; i < size; ++i) {
                hash += seqan3::to_rank(*(query_it)) * fast_pow(k, k - i - 1);
                query_begin++;
            }

            constexpr size_t sigma = seqan3::alphabet_size<alphabet_t>;

            // generate hashes for all kmers that have query as suffix
            std::vector<std::array<int8_t, k>> rank_summands;
            rank_summands.reserve(pow(sigma, k - size));

            std::array<int8_t, k> starter;
            for (size_t i = 0; i < size; ++i)
            {
                if (i >= k - size)
                    starter[i] = seqan3::to_rank(*query_it);
                else
                    starter[i] = -1;

                query_it++;
            }
            rank_summands.push_back(starter);

            std::array<int8_t, sigma> current_n_appended;

            // generate array that holds ranks of chars of generated kmer
            for (int8_t i = k - size - 1; i >= 0; --i)
            {
                // duplicate for next level
                for (auto& h : rank_summands)
                {
                    for (size_t j = 0; j < sigma - 1; ++j)
                    {
                        rank_summands.push_back(h);
                    }
                }

                current_n_appended.fill(0);
                // append in front
                for (auto& h : rank_summands)
                {
                    h[i] = current_n_appended[h.at(i+1)];
                    current_n_appended[h.at(i+1)] += 1;
                }
            }

            // get positions for every kmer with suffix
            std::vector<const std::vector<position_t>*> output;
            for (auto& h : rank_summands)
            {
                size_t current_hash = 0;
                for (size_t i = 0; i < h.size(); ++i)
                    current_hash += h[i] * fast_pow(k, k - i - 1);

                output.emplace_back(_data.find(current_hash));
            }

            // check if query is at the very beginning
            for (size_t i = 0; i < size; ++i)
                if (*backup_start_it != _first_kmer.at(i))
                    return output;

            output.push_back(&_zero);
            return output;
        }

    protected:
        template<std::ranges::range text_t>
        void create(text_t && text)
        {
            // needed for sub_k edge case
            _first_kmer = std::vector<alphabet_t>(text.begin(), text.begin() + k);

            // build map
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

    public:

        template<std::ranges::range text_t>
        minimum_kmer_index(text_t&& text)
        {
            create(std::forward<text_t>(text));
        }

        const std::vector<position_t>& search(std::vector<alphabet_t>& query)
        {
            // search length k directly
            if (query.size() == k)
            {
                const auto& it = _data.find(hash(query.begin()));

                if (it == _data.end())
                    return std::move(std::vector<position_t>{});
                else
                    return it->second;
            }
            else if (query.size() > k) {
                // search first n*k parts
                size_t rest_n = query.size() % k;
                std::vector<std::vector<position_t>*> positions{};
                positions.reserve(query.size() / k + 1);

                for (auto it = query.begin(); it + k != query.end() - rest_n; it += k)
                    positions.push_back(&_data.find(hash(it))->second);

                // search last m < k part
                if (rest_n != 0)
                    for (auto* r : search_subk(query.end() - rest_n, rest_n))
                        positions.push_back(r);

                std::vector<position_t> confirmed_positions{};
                confirmed_positions.reserve(positions.at(0)->size());

                for (auto start_pos : *positions.at(0)) {
                    position_t previous_pos = start_pos;

                    for (size_t i = 1; i <= positions.size(); ++i) {
                        if (i == positions.size()) {
                            confirmed_positions.push_back(start_pos);
                            break;
                        }

                        auto* current = positions.at(i);
                        if (std::find(current->begin(), current->end(), previous_pos + k) != current->end())
                            previous_pos += k;
                        else
                            break;
                    }
                }

                return std::move(confirmed_positions);
            }
            else
            {
                return search_sub_k(query.begin(), query.size());
            }

        }




};