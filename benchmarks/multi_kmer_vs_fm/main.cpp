//
// Created by clem on 6/20/20.
// Copyright (c) 2020 Clemens Cords. All rights reserved.
//
#include <kmer_index.hpp>
#include <benchmarks/input_generator.hpp>
#include <fast_pow.hpp>
#include <benchmarks/cleanup_csv.hpp>

#include <benchmark/benchmark.h>

#include <seqan3/alphabet/nucleotide/dna4.hpp>
#include <seqan3/search/fm_index/fm_index.hpp>
#include <seqan3/search/fm_index/bi_fm_index.hpp>
#include <seqan3/search/search.hpp>

#include <numeric>
#include <random>
#include <memory>

size_t n_benchmarks_registered = 0;
size_t seed = 200;
size_t text_length = 100*1e6;

// [SEQ] exact fm search
template<seqan3::alphabet alphabet_t, typename index_t>
static void fm_search(benchmark::State& state, size_t query_length, size_t text_length, const index_t* index)
{
    state.counters["seed"] = seed;

    input_generator<alphabet_t> input(seed);

    size_t i = 0;
    for (auto _ : state)
    {
        benchmark::DoNotOptimize(search(input.generate_sequence(query_length), *index));
    }

    state.counters["text_length"] = text_length;
    state.counters["query_length"] = query_length;
    state.counters["alphabet_size"] = seqan3::alphabet_size<alphabet_t>;

    seed++;
}

template<seqan3::alphabet alphabet_t, size_t... ks>
static void multi_kmer_search(
        benchmark::State& state,
        size_t query_length,
        size_t text_length,
        const kmer::kmer_index<alphabet_t, uint32_t, ks...>* index)
{
    state.counters["seed"] = seed;
    input_generator<alphabet_t> input(seed);

    bool error_occurred = false;
    try
    {
        size_t i = 0;
        for (auto _ : state)
        {
            benchmark::DoNotOptimize(index->search(input.generate_sequence(query_length)));
        }
    }
    catch (...)
    {
        error_occurred = true;
    }

    state.counters["text_length"] = text_length;
    state.counters["query_length"] = query_length;
    state.counters["alphabet_size"] = seqan3::alphabet_size<alphabet_t>;

    std::vector<size_t> _all_ks = {ks...};
    size_t j = 0;
    for (auto k : _all_ks)
    {
        state.counters["k_" + std::to_string(j)] = k;
        j++;
    }
    state.counters["valid"] = not error_occurred;

    seed++;
}
//nohup ./MULTI_KMER_VS_FM_BENCHMARK --benchmark_format=console --benchmark_counters_tabular=true --benchmark_out=/srv/public/clemenscords/multi_vs_fm/raw.csv --benchmark_out_format=csv --benchmark_repetitions=3 --benchmark_report_aggregates_only=false
int main(int argc, char** argv)
{
    size_t text_length = 1e5;
    auto input = input_generator<seqan3::dna4>(seed);
    auto text = input.generate_sequence(text_length);

    const auto multi_kmer = kmer::make_kmer_index<5, 7, 9, 10, 11, 13, 15, 17, 19, 21, 23, 25, 27, 29, 31>(text, std::thread::hardware_concurrency());
    auto fm = seqan3::fm_index(text);

    text.clear(); // clear memory

    for (size_t query_length = 5; query_length <= 1000; ++query_length)
    {
        benchmark::RegisterBenchmark("multi_kmer", &multi_kmer_search<seqan3::dna4, 5, 7, 9, 10, 11, 13, 15, 17, 19, 21, 23, 25, 27, 29, 31>, query_length, text_length, &multi_kmer);
        benchmark::RegisterBenchmark("fm", &fm_search<seqan3::dna4, decltype(fm)>, query_length, text_length, &fm);
    }

    // start
    benchmark::Initialize(&argc, argv);
    benchmark::RunSpecifiedBenchmarks();

    cleanup_csv("/home/clem/Workspace/kmer_index/source/benchmarks/multi_kmer_vs_fm/new_raw.csv");
}