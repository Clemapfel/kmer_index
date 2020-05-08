// Copyright (c) 2020 Clemens Cords. All rights reserved.

#pragma once

#include <vector>
#include <type_traits>
#include <sstream>

#include <benchmark/benchmark.h>

#include <seqan3/alphabet/concept.hpp>
#include <seqan3/alphabet/aminoacid/all.hpp>
#include <seqan3/alphabet/nucleotide/all.hpp>
#include <seqan3/search/search.hpp>
#include <seqan3/search/fm_index/fm_index.hpp>
#include <seqan3/std/type_traits>

#include "input_generator.hpp"
#include "kmer_index.hpp"

// benchmark arguments struct for readability
template<seqan3::alphabet alphabet_t>
struct benchmark_arguments
{
    private:
        input_generator<alphabet_t> _generator;

        const size_t _query_size;               // size of every query
        const size_t _n_queries;                // number of queries generated
        const size_t _text_size;                // length of text

    public:
        // ctor
        benchmark_arguments(size_t query_size, size_t n_queries, size_t text_size)
                : _query_size(query_size), _n_queries(n_queries), _text_size(text_size)
        {
        }

        // generate queries and text pseudo-randomly
        void generate_queries_and_text(
                std::vector<std::vector<alphabet_t>>* queries,  // out
                std::vector<alphabet_t>* text,                  // out
                bool reset_state = true)                       // if true, returns same queries and text each call
        {
            if (reset_state)
                _generator.reset_state();

            queries->clear();
            for (auto q :_generator.generate_queries(_n_queries, _query_size))
                queries->push_back(q);

            text->clear();
            for (auto c : _generator.generate_text(_text_size, *queries))
                text->push_back(c);
        }

        // add custom counters to keep track of benchmark arguments
        void add_kmer_counters_to(benchmark::State& state, size_t k, bool used_hashtable, size_t n_threads = 1) const
        {
            state.counters["alphabet_size"] = seqan3::alphabet_size<alphabet_t>;
            state.counters["text_size"] = _text_size;
            state.counters["k"] = k;
            state.counters["query_size"] = _query_size;
            state.counters["n_queries"] = _n_queries;
            state.counters["used_hashtable"] = used_hashtable;
            state.counters["paralell"] = n_threads;
        }

        void add_fm_counters_to(benchmark::State& state)
        {
            state.counters["alphabet_size"] = seqan3::alphabet_size<alphabet_t>;
            state.counters["text_size"] = _text_size;
            state.counters["query_size"] = _query_size;
            state.counters["n_queries"] = _n_queries;
        }
};

// #####################################################################################################################

// [SEQ] kmer construction
template<seqan3::alphabet alphabet_t, bool use_da, size_t k>
static void kmer_construction(benchmark::State& state, benchmark_arguments<alphabet_t> input)
{
    std::vector<std::vector<alphabet_t>> queries;
    std::vector<alphabet_t> text;
    input.generate_queries_and_text(&queries, &text, true);

    for (auto _ : state)
        benchmark::DoNotOptimize(debug::make_kmer_index<k>(text, use_da, 1));

    // log memory used
    auto index = debug::make_kmer_index<k>(text, use_da, 1);
    state.counters["memory_used(mb)"] = sizeof(index) / 1e6;

    input.add_counters_to(state, k, use_da);
}

// [PAR] multi kmer construction
template<seqan3::alphabet alphabet_t, bool use_da, size_t... ks>
static void multi_kmer_construction(benchmark::State& state, benchmark_arguments<alphabet_t> input, size_t n_threads)
{
    std::vector<std::vector<alphabet_t>> queries;
    std::vector<alphabet_t> text;
    input.generate_queries_and_text(&queries, &text, true);

    for (auto _ : state)
        benchmark::DoNotOptimize(make_kmer_index<use_da, ks...>(text, n_threads));

    // log memory used
    auto index = make_kmer_index<use_da, ks...>(text, n_threads);
    state.counters["memory_used(mb)"] = sizeof(index) / 1e6;
}

// [SEQ] kmer exact search
template<seqan3::alphabet alphabet_t, bool use_da, size_t k>
static void kmer_search(benchmark::State& state, benchmark_arguments<alphabet_t> input)
{
    // generate text and queries
    std::vector<std::vector<alphabet_t>> queries;
    std::vector<alphabet_t> text;
    input.generate_queries_and_text(&queries, &text, true);

    auto index = debug::make_kmer_index<k>(text, use_da, 1);

    // cycle through queries
    unsigned long long i = 0;
    for (auto _ : state)
    {
        benchmark::DoNotOptimize(index.search(queries.at(i)));
        i = i < queries.size()-1 ? i + 1 : 0;
    }

    input.add_kmer_counters_to(state, k, use_da);
    state.counters["memory_used"] = sizeof(index);
}

// [SEQ] fm exact search
template<seqan3::alphabet alphabet_t>
static void fm_search(benchmark::State& state, benchmark_arguments<alphabet_t> input)
{
    std::vector<std::vector<alphabet_t>> queries;
    std::vector<alphabet_t> text;
    input.generate_queries_and_text(&queries, &text, true);

    seqan3::fm_index index{text};

    unsigned long long i = 0;
    for (auto _ : state)
    {
        benchmark::DoNotOptimize(search(queries.at(i), index));
        i = i < queries.size()-1 ? i + 1 : 0;
    }

    input.add_fm_counters_to(state);
    state.counters["memory_used"] = sizeof(index);
}
