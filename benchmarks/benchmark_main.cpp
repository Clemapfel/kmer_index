// Copyright (c) 2020 Clemens Cords. All rights reserved.

#include "kmer_index.hpp"
#include "benchmarks/input_generator.hpp"
#include "benchmarks/cleanup_csv.hpp"
#include <benchmark/benchmark.h>
#include <seqan3/alphabet/nucleotide/dna4.hpp>
#include <seqan3/search/fm_index/fm_index.hpp>
#include <seqan3/search/search.hpp>

constexpr size_t SEED = 1234;

/* EXECUTABLE ARGUMENTS:
--benchmark_format=console
--benchmark_counters_tabular=true
--benchmark_out=/home/clem/Documents/Workspace/kmer_index/source/benchmarks/fm_kmer_exact_vs_text_length/raw.csv
--benchmark_out_format=csv
--benchmark_repetitions=10
--benchmark_report_aggregates_only=true
 */

template<size_t k>
static void single_kmer_search(benchmark::State& state, size_t text_length, size_t query_length)
{
    auto input = input_generator<seqan3::dna4>(rand());
    auto text = input.generate_sequence(text_length);
    auto queries = input.generate_queries(100000, query_length);

    auto index = kmer::detail::kmer_index_element<seqan3::dna4, k, uint32_t>(text);

    size_t i = 0;
    for (auto _ : state)
    {
        benchmark::DoNotOptimize(index.search(queries.at(i)));
        i = i < queries.size()-1 ? i + 1 : 0;
    }

    state.counters["k"] = k;
    state.counters["text_length"] = text_length;
    state.counters["query_length"] = query_length;
}

template<size_t... ks>
static void multi_kmer_search(benchmark::State& state, size_t text_length, size_t query_length)
{
    auto input = input_generator<seqan3::dna4>(rand());
    auto text = input.generate_sequence(text_length);
    auto queries = input.generate_queries(100000, query_length);

    auto index = kmer_index<seqan3::dna4, uint32_t, ks...>(text);

    size_t i = 0;
    for (auto _ : state)
    {
        benchmark::DoNotOptimize(index.search(queries.at(i)));
        i = i < queries.size()-1 ? i + 1 : 0;
    }

    auto _all_ks = std::vector<size_t>{ks...};

    //state.counters["min_k"] = k;
    state.counters["text_length"] = text_length;
    state.counters["query_length"] = query_length;
}


static void fm_search(benchmark::State& state, size_t text_length, size_t query_length)
{
    auto input = input_generator<seqan3::dna4>(rand());
    auto text = input.generate_sequence(text_length);
    auto queries = input.generate_queries(100000, query_length);

    auto index = seqan3::fm_index(text);

    size_t i = 0;
    for (auto _ : state)
    {
        benchmark::DoNotOptimize(seqan3::search(queries.at(i), index));
        i = i < queries.size()-1 ? i + 1 : 0;
    }

    state.counters["text_length"] = text_length;
    state.counters["query_length"] = query_length;
}

template<size_t k>
void register_benchmarks(size_t text_length, size_t query_length)
{
    auto input = input_generator<seqan3::dna4>(SEED);
    auto text = input.generate_sequence(text_length);
    auto queries = input.generate_queries(query_length, 100000);

    benchmark::RegisterBenchmark("single", &single_kmer_search<k>, text_length, query_length);
    benchmark::RegisterBenchmark("multi", &multi_kmer_search<k>, text_length, query_length);
}

constexpr size_t k = 6;
constexpr size_t k_2 = k - 1;
constexpr size_t k_4 = k + 1;

int main(int argc, char** argv)
{
    //register_benchmarks<6>(100000, k);
    //register_benchmarks<6>(100000, k-2);

    for (size_t i = k-2; i < 3*k; i += 1)
        register_benchmarks<k>(100000, i);

    benchmark::Initialize(&argc, argv);
    benchmark::RunSpecifiedBenchmarks();

    cleanup_csv("/home/clem/Documents/Workspace/kmer_index/benchmarks/multi_vs_single/raw.csv");
}