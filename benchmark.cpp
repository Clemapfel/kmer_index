//
// Created by Clemens Cords on 2/7/20.
// Copyright (c) 2020 Clemens Cords. All rights reserved.
//
#include <vector>

#include <benchmark/benchmark.h>

#include <seqan3/alphabet/concept.hpp>
#include <seqan3/search/algorithm/search.hpp>
#include <seqan3/search/fm_index/fm_index.hpp>

#include "input_generator.hpp"
#include "kmer_index.hpp"

// register programmatically with registerbenchmark: https://github.com/google/benchmark#output-files

/*
should_use_DA
k
alphabet
n_queries
query_size
query_size_distribution: half non n*k, half nk for average case
n hits per query
% of queries with 0 hits

constexpr size_t
 */

// input config in a single struct for convenience
struct input_generator_config
{
    input_generator_config(
        size_t n_queries,                   // total number of queries
        std::vector<size_t> query_lengths,  // possible query lengths
        size_t text_length,                 // length of text
        float no_hit_queries_ratio = 0)     // percentage of queries that will have no occurence in text
            : _n_queries(n_queries),
              _query_lengths(query_lengths),
              _text_length(text_length),
              _no_hit_queries_ratio(no_hit_queries_ratio)
    {
    }

    const std::vector<size_t> _query_lengths;
    const size_t _n_queries;
    const size_t _text_length;
    const float _no_hit_queries_ratio;

    std::string name()
    {
        return "todo";
    }
};

// fill queries and vector based on config
template<seqan3::alphabet alphabet_t>
void generate_queries_and_text(
    input_generator_config config,                  // in
    std::vector<std::vector<alphabet_t>>* queries,  // out
    std::vector<alphabet_t>* text)                  // out
{
    input_generator<alphabet_t>::reset_state();

    queries->clear();
    for (auto q : input_generator<alphabet_t>::generate_queries(config._n_queries, config._query_lengths))
        queries->push_back(q);

    text->clear();
    for (auto c : input_generator<alphabet_t>::generate_text(config._text_length, *queries))
        text->push_back(c);
}

// benchmark: kmer_index exact search time
template<seqan3::alphabet alphabet_t, bool use_da, size_t k>
static void kmer_search(benchmark::State& state)
{
    auto input = input_generator_config(10, std::vector<size_t>{5, 6, 7}, 1000, 0.f);

    // generate text and queries
    std::vector<std::vector<alphabet_t>> queries;
    std::vector<alphabet_t> text;
    generate_queries_and_text(input, &queries, &text);

    kmer_index<alphabet_t, k, uint64_t, uint32_t, use_da> index{text};

    // cycle through queries
    unsigned long long i = 0;
    for (auto _ : state)
    {
        index.search(queries.at(i));
        i = i < queries.size()-1 ? i + 1 : 0;
    }
}

// benchmark: fm_index exact search time
template<seqan3::alphabet alphabet_t>
static void fm_search(benchmark::State& state, input_generator_config input)
{
    // generate text and queries
    std::vector<std::vector<alphabet_t>> queries;
    std::vector<alphabet_t> text;
    generate_queries_and_text(input, &queries, &text);

    seqan3::fm_index index{text};

    // cycle through queries
    unsigned long long i = 0;
    for (auto _ : state)
    {
        seqan3::search(index, queries.at(i));
        i = i < queries.size() ? i + 1 : 0;
    }
}

// wrapper that programmatically registers all
template<seqan3::alphabet alphabet_t, bool use_da, size_t... ks>
struct register_benchmarks
{
    static void call(std::vector<input_generator_config> configs)
    {
        for (auto config : configs)
        {
            (benchmark::RegisterBenchmark(config.name().c_str(), &kmer_search<alphabet_t, use_da, ks>), ...);
        }

        benchmark::Initialize(0, nullptr);
        benchmark::RunSpecifiedBenchmarks();
    }

        const std::vector<size_t> _query_lengths;
    const size_t _n_queries;
    const size_t _text_length;
    const float _no_hit_queries_ratio;
};

int main()
{
    auto cfg = input_generator_config(10, std::vector<size_t>{5, 6, 7}, 1000, 0.f);
    register_benchmarks<seqan3::dna4, false, 4, 5, 6, 7, 8>::call({cfg});

}



