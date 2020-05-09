// Copyright (c) 2020 Clemens Cords. All rights reserved.

#include <vector>
#include <type_traits>
#include <sstream>

#include <benchmark/benchmark.h>

#include <seqan3/alphabet/concept.hpp>
#include <seqan3/alphabet/aminoacid/all.hpp>
#include <seqan3/alphabet/nucleotide/all.hpp>
#include <seqan3/search/search.hpp>
#include <seqan3/search/fm_index/fm_index.hpp>

#include "input_generator.hpp"
#include "kmer_index.hpp"
#include "cleanup_csv.hpp"

struct benchmark_arguments
{
    size_t text_length;
    size_t query_length;
    size_t n_queries;

    void add_counters_to(benchmark::State& state)
    {
        state.counters["text_length"] = text_length;
        state.counters["query_length"] = query_length;
        state.counters["n_queries"] = n_queries;
    }
};

constexpr size_t SEED = 9999;

// [SEQ] kmer exact search
template<seqan3::alphabet alphabet_t, size_t k>
static void kmer_search(benchmark::State& state, benchmark_arguments args)
{
    auto input = input_generator<alphabet_t>(SEED);

    // generate text
    std::vector<alphabet_t> text = input.generate_sequence(args.text_length);

    // generate queries
    std::vector<std::vector<alphabet_t>> queries;
    for (size_t i = 0; i < args.n_queries; ++i)
        queries.push_back(input.generate_sequence(args.query_length));

    auto index = minimum_kmer_index<alphabet_t, k>(text);

    // cycle through queries
    unsigned long long i = 0;
    for (auto _ : state)
    {
        benchmark::DoNotOptimize(index.search(queries.at(i)));
        i = i < queries.size()-1 ? i + 1 : 0;
    }

    state.counters["k"] = k;
    args.add_counters_to(state);
}

// [SEQ] fm exact search
template<seqan3::alphabet alphabet_t>
static void fm_search(benchmark::State& state, benchmark_arguments args)
{
    auto input = input_generator<alphabet_t>(SEED);

    using namespace seqan3;

    // generate text
    std::vector<alphabet_t> text = input.generate_sequence(args.text_length);

    // generate queries
    std::vector<std::vector<alphabet_t>> queries;
    for (size_t i = 0; i < args.n_queries; ++i)
        queries.push_back(input.generate_sequence(args.query_length));

    seqan3::fm_index index{text};

    unsigned long long i = 0;
    for (auto _ : state)
    {
        benchmark::DoNotOptimize(search(queries.at(i), index));
        i = i < queries.size()-1 ? i + 1 : 0;
    }

    args.add_counters_to(state);
}

inline size_t N_BENCHMARKS = 0;

// register all benchmarks
template<seqan3::alphabet alphabet_t, size_t k>
void register_benchmarks(size_t text_size)
{
    for (size_t i = 1; i < 10; ++i)
    {
        auto args = benchmark_arguments{text_size, i * k, 100000};

        benchmark::RegisterBenchmark("kmer_search", &kmer_search<alphabet_t, k>, args);
        benchmark::RegisterBenchmark("fm_search", &fm_search<alphabet_t>, args);
        N_BENCHMARKS++;
    }
}

template<size_t... ks>
void register_all(size_t text_length, int argc, char** argv)
{
    (register_benchmarks<seqan3::dna4, ks>(text_length), ...);
}

/* EXECUTABLE ARGUMENTS:

 *
--benchmark_format=console
--benchmark_counters_tabular=true
--benchmark_out=/home/clem/Documents/Workspace/kmer_index/source/benchmarks/minimum_raw.csv
--benchmark_out_format=csv
--benchmark_repetitions=10
--benchmark_report_aggregates_only=true
 */

int main(int argc, char** argv)
{
    register_all<4, 5, 6, 7, 8>(500000, argc, argv);

    benchmark::Initialize(&argc, argv);
    seqan3::debug_stream << std::to_string(N_BENCHMARKS) + " benchmarks currently registered\n";
    benchmark::RunSpecifiedBenchmarks();

    //cleanup_csv("/home/clem/Documents/Workspace/kmer_index/source/benchmarks/minimum_raw.csv");
}