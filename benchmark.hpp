//
// Created by Clemens Cords on 3/3/20.
// Copyright (c) 2020 Clemens Cords. All rights reserved.
//

#pragma once

#include <vector>
#include <type_traits>
#include <sstream>

#include <benchmark/benchmark.h>

#include <seqan3/alphabet/concept.hpp>
#include <seqan3/alphabet/aminoacid/all.hpp>
#include <seqan3/alphabet/nucleotide/all.hpp>
#include <seqan3/search/algorithm/search.hpp>
#include <seqan3/search/fm_index/fm_index.hpp>
#include <seqan3/std/type_traits>

#include "input_generator.hpp"
#include "kmer_index.hpp"


// benchmark arguments struct for readability
struct benchmark_arguments
{
    const size_t _query_size;               // size of every query
    const size_t _n_queries;                // number of queries generated
    const size_t _text_size;                // length of text

    // ctor
    benchmark_arguments(size_t query_size, size_t n_queries, size_t text_size)
        : _query_size(query_size), _n_queries(n_queries), _text_size(text_size)
    {
    }

    // generate queries and text pseudo-randomly
    template<seqan3::alphabet alphabet_t>
    void generate_queries_and_text(
        std::vector<std::vector<alphabet_t>>* queries,  // out
        std::vector<alphabet_t>* text,                  // out
        bool reset_state = true)                       // if true, returns same queries and text each call
    {
        if (reset_state)
            input_generator<alphabet_t>::reset_state();

        queries->clear();
        for (auto q : input_generator<alphabet_t>::generate_queries(_n_queries, _query_size))
            queries->push_back(q);

        text->clear();
        for (auto c : input_generator<alphabet_t>::generate_text(_text_size, *queries))
            text->push_back(c);
    }

    // add custom counters to keep track of benchmark arguments
    template<seqan3::alphabet alphabet_t, bool use_hashtable, int k>
    void add_counters_to(benchmark::State& state) const
    {
        state.counters["alphabet_size"] = seqan3::alphabet_size<alphabet_t>;
        state.counters["text_size"] = _text_size;
        state.counters["k"] = k;
        state.counters["query_size"] = _query_size;
        state.counters["n_queries"] = _n_queries;
        state.counters["used_hashtable"] = use_hashtable;
    }
};

// #####################################################################################################################

// benchmark: kmer index exact search
template<seqan3::alphabet alphabet_t, bool use_da, size_t k>
static void kmer_search(benchmark::State& state, benchmark_arguments input)
{
    // generate text and queries
    std::vector<std::vector<alphabet_t>> queries;
    std::vector<alphabet_t> text;
    input.generate_queries_and_text(&queries, &text, true);

    auto index = make_kmer_index<use_da, k>(text);

    // add custom vars to output
    input.add_counters_to<alphabet_t, use_da, k>(state);
    state.counters["memory_used(mb)"] = sizeof(index) / 1e6;

    // cycle through queries
    unsigned long long i = 0;
    for (auto _ : state)
    {
        benchmark::DoNotOptimize(index.search(queries.at(i)));
        i = i < queries.size()-1 ? i + 1 : 0;
    }
}

/*
// benchmark:: multi kmer index exact search
template<seqan3::alphabet alphabet_t, bool use_da, size_t... ks>
static void multi_kmer_search(benchmark::State& state, benchmark_arguments input)
{
    std::vector<std::vector<alphabet_t>> queries;
    std::vector<alphabet_t> text;
    input.generate_queries_and_text(&queries, &text, true);

    multi_kmer_index<alphabet_t, k, uint64_t, uint32_t, use_da, ks...> index{text};

    unsigned long long i = 0;
    for (auto _ : state)
    {
        benchmark::DoNotOptimize(index.search(queries.at(i)));
        i = i < queries.size()-1 ? i + 1 : 0;
    }
}
*/

// benchmark: fm index exact search
template<seqan3::alphabet alphabet_t>
static void fm_search(benchmark::State& state, benchmark_arguments input)
{
    std::vector<std::vector<alphabet_t>> queries;
    std::vector<alphabet_t> text;
    input.generate_queries_and_text(&queries, &text, true);

    seqan3::fm_index index{text};

    input.add_counters_to<alphabet_t, false, -1>(state);
    state.counters["memory_used(mb)"] = sizeof(index) / 1e6;

    unsigned long long i = 0;
    for (auto _ : state)
    {
        benchmark::DoNotOptimize(search(queries.at(i), index));
        i = i < queries.size()-1 ? i + 1 : 0;
    }
}

// benchmark: kmer_index construction
template<seqan3::alphabet alphabet_t, bool use_da, size_t k>
static void kmer_construction(benchmark::State& state, benchmark_arguments input)
{
    std::vector<std::vector<alphabet_t>> queries;
    std::vector<alphabet_t> text;
    input.generate_queries_and_text<seqan3::dna4>(&queries, &text, true);

    for (auto _ : state)
        benchmark::DoNotOptimize(make_kmer_index<use_da, k>(text));

    // log memory used
    auto index = make_kmer_index<use_da, k>(text);
    state.counters["memory_used(mb)"] = sizeof(index) / 1e6;
}

/*
// benchmark: multi_kmer construction
template<seqan3::alphabet alphabet_t, bool use_da, size_t... ks>
static void multi_kmer_search(benchmark::State& state, benchmark_arguments input)
{
    std::vector<std::vector<alphabet_t>> queries;
    std::vector<alphabet_t> text;-
    input.generate_queries_and_text(&queries, &text, true);

    for (auto _ : state)
        benchmark::DoNotOptimize(multi_kmer_index<alphabet_t, k, uint64_t, uint32_t, use_da, ks...> index{text});
}
 */

// benchmark: fm_index construction
template<seqan3::alphabet alphabet_t>
static void fm_construction(benchmark::State& state, benchmark_arguments input)
{
    std::vector<std::vector<alphabet_t>> queries;
    std::vector<alphabet_t> text;
    input.generate_queries_and_text<seqan3::dna4>(&queries, &text, true);

    for (auto _ : state)
        benchmark::DoNotOptimize(seqan3::fm_index{text});

    // log memory used
    seqan3::fm_index index{text};
    input.add_counters_to<alphabet_t, false, -1>(state);
    state.counters["memory_used(mb)"] = sizeof(index) / 1e6;
}


//#####################################################################################################################

// function that can be expanded in fold expression
template<seqan3::alphabet alphabet_t, bool use_da, size_t k>
void register_kmer_benchmarks(benchmark_arguments config)
{
    benchmark::RegisterBenchmark("kmer_exact_search", &kmer_search<alphabet_t, use_da, k>, config);
    benchmark::RegisterBenchmark("kmer_construction", &kmer_construction<alphabet_t, use_da, k>, config);;
}

// wrapper that programmatically registers benchmarks for all permutations of template params and args
template<seqan3::alphabet alphabet_t, bool use_hashtable, size_t... ks>
void register_all_benchmarks(std::vector<size_t> text_sizes, int query_size_offset = 0)
{
    assert(k - query_size > 0);

    std::map<size_t, std::vector<benchmark_arguments>> _configs;

    for (auto k : std::vector{ks...})
    {
        size_t query_size = k + query_size_offset;

        for (auto text_size : text_sizes)
        {
            size_t n_queries = 0.1 * text_size > 100 ? 0.1 * text_size : 100;

            if (_configs.find(k) == _configs.end())
                _configs.insert(std::make_pair(k, std::vector<benchmark_arguments>{}));

            _configs[k].emplace_back(query_size, n_queries, text_size);
        }
    }

    for (auto pair : _configs)
    {
        for (auto& config : pair.second)
        {
            (register_kmer_benchmarks<alphabet_t, use_hashtable, ks>(config), ...);
            benchmark::RegisterBenchmark("fm_search", &fm_search<alphabet_t>, config);
            benchmark::RegisterBenchmark("fm_construction", &fm_construction<alphabet_t>, config);
        }
    }
}

// cleanup resulting csv for later parsing
void cleanup_csv(std::string path)
{
    std::ifstream file_in;
    file_in.open(path);

    if (file_in.fail())
        std::cerr << "Error opening file " + path + "\naborting...";

    // name outfile
    auto it = path.find("_raw");
    std::string out_path = path.erase(it, 4);

    std::ofstream file_out;
    file_out.open(out_path);

    if (file_out.fail())
        std::cerr << "Error opening file " + out_path + "\naborting...";

    // discard and reformat header
    std::string line;
    bool past_header = false;

    while (std::getline(file_in, line))
    {
        if (not past_header)
        {
            // csv header line
            if (line.substr(0, 4) == "name")
            {
                line.erase(std::remove(line.begin(), line.end(), '\"'), line.end());
                file_out << line + "\n";
                past_header = true;
            }
            // skip gbenchmark header
            else
                continue;
        }
        else
            file_out << line + "\n";
    }

    std::cout << "[DEBUG] csv written to " + out_path + "\n";

    file_in.close();
    file_out.close();
}
