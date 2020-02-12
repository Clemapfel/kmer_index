//
// Created by Clemens Cords on 2/7/20.
// Copyright (c) 2020 Clemens Cords. All rights reserved.
//
#include <vector>
#include <sstream>

#include <benchmark/benchmark.h>

#include <seqan3/alphabet/concept.hpp>
#include <seqan3/search/algorithm/search.hpp>
#include <seqan3/search/fm_index/fm_index.hpp>

#include "input_generator.hpp"
#include "kmer_index.hpp"

// benchmark arguments struct for readability
struct benchmark_arguments
{
    const size_t _query_size;               // size of every query
    const size_t _n_queries;                // number of queries generated
    const size_t _text_size;                // length of text
    const float _no_hit_queries_ratio;      // text containts ratio * n_queries of the queries, 0 hits for rest

    // ctor
    benchmark_arguments(
        size_t query_size,
        size_t n_queries,
        size_t text_size,
        float no_hit_queries_ratio = 0)
        : _query_size(query_size),
          _n_queries(n_queries),
          _text_size(text_size),
          _no_hit_queries_ratio(no_hit_queries_ratio)
    {
    }

    // generate queries and text pseudo-randomly
    template<seqan3::alphabet alphabet_t>
    void generate_queries_and_text(
        std::vector<std::vector<alphabet_t>* queries, // out
        std::vector<alphabet_t> text,                 // out
        bool reset_state = true)                      // if true, returns same queries and text each call
    {
        if (reset_state)
            input_generator<alphabet_t>::reset_state();

        if (queries != nullptr)
        {
            queries->clear();
            for (auto q : input_generator<alphabet_t>::generate_queries(config._n_queries, config._query_lengths))
                queries->push_back(q);
        }

        if (text != nullptr)
        {
            text->clear();
            for (auto c : input_generator<alphabet_t>::generate_text(config._text_length, *queries))
                text->push_back(c);
        }
    }

     // generate csv header
    static std::string get_header()
    {
        return "alphabet,text_size,k,query_size,n_queries,no_hit_query_ratio,uses_direct_addressing,";
    }

    // generate benchmark name which forms first few columns of csv
    std::string get_name(size_t k, std::string alphabet_id, bool use_da)
    {
        return
           "\"" + alphabet_id + "\","
           + std::to_string(_text_size) + ","
           + std::to_string(k) + ","
           + std::to_string(_query_size) + ","
           + std::to_string(_n_queries )+ ","
           + std::to_string(_no_hit_queries_ratio) + ","
           + std::to_string(use_da) + ",";
    }
};

// benchmark: kmer index exact search
template<seqan3::alphabet alphabet_t, bool use_da, size_t k>
static void kmer_search(benchmark::State& state, benchmark_arguments input)
{
    // generate text and queries
    std::vector<std::vector<alphabet_t>> queries;
    std::vector<alphabet_t> text;
    input.generate_queries_and_text(&queries, &text, true);

    kmer_index<alphabet_t, k, uint64_t, uint32_t, use_da> index{text};

    // cycle through queries
    unsigned long long i = 0;
    for (auto _ : state)
    {
        index.search(queries.at(i));
        i = i < queries.size()-1 ? i + 1 : 0;
    }
}

// benchmark: fm index exact search
template<seqan3::alphabet alphabet_t, bool use_da, size_t k>
static void fm_search(benchmark::State& state, benchmark_arguments input)
{
    std::vector<std::vector<alphabet_t>> queries;
    std::vector<alphabet_t> text;
    input.generate_queries_and_text(&queries, &text, true);

    fm_index index{text};

    unsigned long long i = 0;
    for (auto _ : state)
    {
        search(index, queries.at(i));
        i = i < queries.size()-1 ? i + 1 : 0;
    }
}

// benchmark: kmer_index construction
template<seqan3::alphabet alphabet_t, bool use_da, size_t k>
static void kmer_construction(benchmark::State& state, benchmark_arguments input)
{
    std::vector<alphabet_t> text;
    input.generate_queries_and_text(nullptr, &text, true);

    for (auto _ : state)
        benchmark::DoNotOptimize(kmer_index<alphabet_t, k, uint64_t, uint32_t, use_da> index{text});

    // log memory used
    state["memory_used"] = kmer_index<alphabet_t, k, uint64_t, uint32_t, use_da>::estimate
}

// main
int main(int argc, char** argv) {

    auto cfg = benchmark_arguments(100, 10, 1000, 0);
    std::cout << cfg.get_name(5, "dna4", true);
}


