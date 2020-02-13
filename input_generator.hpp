//
// Created by Clemens Cords on 2/7/20.
// Copyright (c) 2020 Clemens Cords. All rights reserved.
//

#pragma once

#include <random>
#include <vector>

#include <seqan3/alphabet/concept.hpp>
#include <seqan3/alphabet/hash.hpp>
#include <seqan3/range/views/kmer_hash.hpp>

// input generator providing high quality pseudo-random sequences
template<seqan3::alphabet alphabet_t, size_t seed = 1234>
class input_generator
{
    private:
        // all function calls use this engine and advance it's state
        inline static std::mt19937 _engine{seed};

        static uint64_t hash(std::vector<alphabet_t> query)
        {
            if (query.size() == 0)
                return 0;

            auto shape = seqan3::shape{seqan3::ungapped{static_cast<uint8_t>(query.size())}};
            return static_cast<uint64_t>( *(query | seqan3::views::kmer_hash(shape)).begin());
        }

    public:
        // resets rng state
        // can be used before generation call to get the same results as when the class was first used
        static void reset_state()
        {
            _engine.seed(seed);
        }

        // generate sequence
        static std::vector<alphabet_t> generate_sequence(size_t length)
        {
            std::uniform_int_distribution<uint8_t> dist(0, seqan3::alphabet_size<alphabet_t> -1);
            std::vector <alphabet_t> sequence;

            for (size_t l = 0; l < length; ++l)
            {
                auto number = dist(_engine);
                sequence.push_back(alphabet_t{}.assign_rank(number));
            }
            return sequence;
        }

        // overload: generate a total of n queries of size length
        static std::vector<std::vector<alphabet_t>> generate_queries(size_t n_queries, size_t length)
        {
            std::vector<std::vector<alphabet_t>> queries;

            for (size_t i = 0; i < n_queries; ++i)
                queries.push_back(generate_sequence(length));

            return queries;
        }

        // generate text that contains queries
        static std::vector<alphabet_t> const generate_text(size_t length, std::vector<std::vector<alphabet_t>> queries, float missing_ratio = 0)
        {
            size_t max_length = 0;
            for (auto q : queries)
                if (q.size() > max_length)
                    max_length = q.size();

            // toss out percentage of queries
            size_t n_tossed_out = missing_ratio * queries.size();

            std::unordered_set<size_t> missing{};

            for (size_t i = 0; i < n_tossed_out; ++i)
            {
                std::uniform_int_distribution<size_t> which_query(0, queries.size() - 1);
                size_t which = which_query(_engine);
                missing.insert(input_generator::hash(queries.at(which)));
                queries.erase(queries.begin() + which);
            }

            // generate
            std::uniform_real_distribution<float> rng_chance(0.01, 0.1);
            auto query_insert_chance = rng_chance(_engine);

            std::bernoulli_distribution bernoulli(query_insert_chance);
            std::uniform_int_distribution<size_t> which_query(0, queries.size() - 1);
            std::uniform_int_distribution<size_t> random_insert_length(0, max_length);

            std::vector<alphabet_t> text{};

            while (text.size() < length)
            {
                if (bernoulli(_engine))
                {
                    for (auto c : queries.at(which_query(_engine)))
                        text.push_back(c);
                }
                else
                {
                    auto to_insert = generate_sequence(random_insert_length(_engine));
                    for (auto &c : to_insert)
                        text.push_back(c);
                }
            }

            text.resize(length);

            return text;
        }
};