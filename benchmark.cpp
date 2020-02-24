//
// Created by Clemens Cords on 2/7/20.
// Copyright (c) 2020 Clemens Cords. All rights reserved.
//
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
        std::vector<std::vector<alphabet_t>>* queries,  // out
        std::vector<alphabet_t>* text,                  // out
        bool reset_state = true)                        // if true, returns same queries and text each call
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
    template<seqan3::alphabet alphabet_t, bool use_da, int k>
    void add_counters_to(benchmark::State& state) const
    {
        state.counters["alphabet_size"] = seqan3::alphabet_size<alphabet_t>;
        state.counters["text_size"] = _text_size;
        state.counters["k"] = k;
        state.counters["query_size"] = _query_size;
        state.counters["n_queries"] = _n_queries;
        state.counters["no_hit_query_ratio"] = _no_hit_queries_ratio;
        state.counters["used_da"] = use_da;
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

    auto index = make_kmer_index<k>(text);

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
        benchmark::DoNotOptimize(make_kmer_index<k, use_da>(text));

    // log memory used
    auto index = make_kmer_index<k>(text);
    input.add_counters_to<alphabet_t, use_da, k>(state);
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
bool register_kmer_benchmarks(benchmark_arguments config)
{
    // skip invalid config automatically
    if (detail::pow_ul(seqan3::alphabet_size<alphabet_t>, k) < UINT64_MAX)
        return false;

    benchmark::RegisterBenchmark("kmer_exact_search", &kmer_search<alphabet_t, use_da, k>, config);
    benchmark::RegisterBenchmark("kmer_construction", &kmer_construction<alphabet_t, use_da, k>, config);

    return true;
}

// wrapper that programmatically registers benchmarks for all permutations of template params and args
template<seqan3::alphabet alphabet_t, bool use_da, size_t... ks>
void register_all_benchmarks(
    std::map<size_t, std::vector<size_t>> query_sizes_per_k,
    std::vector<size_t> text_sizes,
    std::vector<float> no_hit_ratios)
{
    std::map<size_t, std::vector<benchmark_arguments>> _configs;

    for (auto k : std::vector{ks...})
    {
        std::vector<size_t> query_sizes;
        if (query_sizes_per_k.find(k) != query_sizes_per_k.end())
            query_sizes = query_sizes_per_k.at(k);
        else
            query_sizes = {k};

        for (auto query_size : query_sizes)
        {
            for (auto text_size : text_sizes)
            {
                for (auto no_hit_ratio : no_hit_ratios)
                {
                    size_t n_queries = 0.1 * text_size > 100 ? 0.1 * text_size : 100;

                    if (_configs.find(k) == _configs.end())
                        _configs.insert(std::make_pair(k, std::vector<benchmark_arguments>{}));

                    _configs[k].emplace_back(query_size, n_queries, text_size, no_hit_ratio);
                }
            }
        }
    }

    for (auto pair : _configs)
    {
        for (auto& config : pair.second)
        {
            (register_kmer_benchmarks<alphabet_t, use_da, ks>(config, true), ...);
            benchmark::RegisterBenchmark("fm_search", &fm_search<alphabet_t>, config);
            benchmark::RegisterBenchmark("fm_construction", &fm_construction<alphabet_t>, config);
        }
    }
}

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
                file_out << line;
                past_header = true;
            }
            // skip gbenchmark header
            else
                continue;
        }
        else
            file_out << line;
    }

    std::cout << "[DEBUG] csv written to " + out_path + "\n";

    file_in.close();
    file_out.close();
}

// #####################################################################################################################

std::index_sequence<5, 10, 15, 20, 25, 30> all_ks;
std::vector<size_t> text_sizes = {50e3, 100e3, 250e3, 500e3, 1e6};//, 2e6, 5e6, 10e6, 250e6, 3e9};
std::vector<float> no_hit_ratio = {0, 0.25, 0.5, 0.75, 1};

// main
int main(int argc, char** argv)
{
    std::map<size_t, std::vector<size_t>> query_sizes;
    query_sizes[20] = {19, 20};
    query_sizes[25] = {24, 25};

    register_all_benchmarks<seqan3::dna4, false, all_ks>(query_sizes, text_sizes, hit_ratios);
    register_all_benchmarks<seqan3::dna4, true, all_ks>(query_sizes, text_sizes, hit_ratios);

    benchmark::Initialize(&argc, argv);
    benchmark::RunSpecifiedBenchmarks();

    cleanup_csv("../source/benchmark_out_raw.csv");
}


