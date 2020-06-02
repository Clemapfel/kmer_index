// Copyright (c) 2020 Clemens Cords. All rights reserved.

#pragma once

#include <seqan3/alphabet/nucleotide/dna4.hpp>
#include <seqan3/search/fm_index/fm_index.hpp>
#include <seqan3/search/search.hpp>

#include <benchmark/benchmark.h>

#include "kmer_index.hpp"
#include "benchmarks/input_generator.hpp"
#include "benchmarks/cleanup_csv.hpp"


// benchmark input struct that holds arguments used for generating texts
    template<seqan3::alphabet alphabet_t>
    struct benchmark_input
    {
            // CTOR
            benchmark_input(size_t text_length, size_t query_length, size_t seed = rand())
                    : _text_length(text_length), _query_length(query_length), _input(seed)
            {
            }

            // generate text
            std::vector<alphabet_t> generate_text() const
            {
                return _input.generate_sequence(_text_length);
            }

            // queries
            std::vector<std::vector<alphabet_t>> generate_queries() const
            {
                std::vector<std::vector<alphabet_t>> output;
                output.reserve(_query_length);

                for (size_t i = 0; i < _n_queries; ++i)
                    output.push_back(_input.generate_sequence(_query_length));

                return output;
            }

            // add counters to benchmark csv
            template<size_t... ks>
            void add_counters_to(benchmark::State& state)
            {
                /*
                std::vector<size_t> all_ks = {(ks, ...)};
                for (size_t i = 0; i < all_ks.size(); ++i)
                    state.counters["k_" + std::to_string(i)] = all_ks.at(i);

                state.counters["text_length"] = _text_length;
                state.counters["query_length"] = _query_length;
                state.counters["alphabet_size"] = seqan3::alphabet_size<alphabet_t>;*/
            }

            // overload for fm_index which doesn't have k
            void add_counters_to(benchmark::State& state)
            {
                state.counters["text_length"] = _text_length;
                state.counters["query_length"] = _query_length;
                state.counters["alphabet_size"] = seqan3::alphabet_size<alphabet_t>;
            }

        private:
            const size_t _text_length;
            const size_t _query_length;

            mutable input_generator<alphabet_t> _input;

            inline static const size_t _n_queries = 100000;
    };

// ###############################################################################################

// [SEQ] exact single kmer search
        template<seqan3::alphabet alphabet_t, size_t k>
        static void kmer_search(benchmark::State& state, benchmark_input<alphabet_t> input)
        {
            using namespace seqan3;
            auto text = input.generate_text();
            std::vector<std::vector<dna4>> queries = input.generate_queries();

            auto index = kmer::make_kmer_index<k>(text);

            size_t i = 0;
            for (auto _ : state)
            {
                benchmark::DoNotOptimize(index.search(queries.at(i)));
                i = i < queries.size() - 1 ? i + 1 : 0;
            }

            input.add_counters_to(state);
            state.counters["k"] = k;
        }

// [SEQ] exact multi kmer search
    template<seqan3::alphabet alphabet_t, size_t... ks>
    static void multi_kmer_search(benchmark::State& state, benchmark_input<alphabet_t> input)
    {
        auto text = input.generate_sequence();
        auto queries = input.generate_queries();

        auto index = kmer::make_kmer_index<ks...>(text);

        size_t i = 0;
        for (auto _ : state)
        {
            benchmark::DoNotOptimize(index.search(queries.at(i)));
            i = i < queries.size() - 1 ? i + 1 : 0;
        }

        input.add_counters_to<(ks, ...)>(state);
    }


// [SEQ] exact fm index search
    template<seqan3::alphabet alphabet_t>
    static void fm_search(benchmark::State& state, benchmark_input<alphabet_t> input)
    {
        auto text = input.generate_text();
        auto queries = input.generate_queries();

        auto index = seqan3::fm_index(text);

        size_t i = 0;
        for (auto _ : state)
        {
            benchmark::DoNotOptimize(seqan3::search(queries.at(i), index));
            i = i < queries.size() - 1 ? i + 1 : 0;
        }

        input.add_counters_to(state);
    }


