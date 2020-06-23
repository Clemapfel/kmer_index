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

size_t n_benchmarks_registered = 0;
size_t seed = 200;
size_t text_length = 1000000;

// [SEQ] exact fm search
template<seqan3::alphabet alphabet_t>
static void fm_search(benchmark::State& state, size_t query_length, const std::vector<alphabet_t>& text)
{
    state.counters["seed"] = seed;

    input_generator<alphabet_t> input(seed);

    auto index = seqan3::fm_index(text);

    size_t i = 0;
    for (auto _ : state)
    {
        benchmark::DoNotOptimize(search(input.generate_sequence(query_length), index));
    }

    state.counters["text_length"] = text_length;
    state.counters["query_length"] = query_length;
    state.counters["alphabet_size"] = seqan3::alphabet_size<alphabet_t>;

    seed++;
}

template<seqan3::alphabet alphabet_t, size_t k>
static void single_kmer_search(
        benchmark::State& state,
        size_t query_length,
        const kmer::kmer_index<alphabet_t, uint32_t, k>* index)
{
    state.counters["seed"] = seed;
    input_generator<alphabet_t> input(seed);

    bool error_occurred = false;
    try
    {
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
    state.counters["k_0"] = k;
    state.counters["valid"] = not error_occurred;

    seed++;
}

template<seqan3::alphabet alphabet_t, size_t... ks>
static void multi_kmer_search(
        benchmark::State& state,
        size_t query_length,
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

void register_all()
{
    auto input = input_generator<seqan3::dna4>(seed);
    const auto text = input.generate_sequence(text_length);

    std::cout << "starting building\n";

    const auto multi_kmer = kmer::make_kmer_index<10, 11, 13, 15, 17, 19, 21, 23, 25, 27>(text, std::thread::hardware_concurrency());
    const auto single_kmer = kmer::make_kmer_index<15>(text, 1);
    auto fm = seqan3::fm_index(text);

    std::cout << "building finished.\n";

    for (size_t query_length = 7; query_length <= 200; ++query_length)
    {
        benchmark::RegisterBenchmark("multi_kmer", &multi_kmer_search<seqan3::dna4, 10, 11, 13, 15, 17, 19, 21, 23, 25, 27>, query_length, &multi_kmer);
        benchmark::RegisterBenchmark("single_kmer", &single_kmer_search<seqan3::dna4, 15>, query_length, &single_kmer);
        benchmark::RegisterBenchmark("fm", &fm_search<seqan3::dna4>, query_length, text);

        n_benchmarks_registered += 3;
    }

    std::cout << std::to_string(n_benchmarks_registered) << " benchmarks registered.\n";
}

//nohup ./MULTI_KMER_BENCHMARK --benchmark_format=console --benchmark_counters_tabular=true --benchmark_out=/srv/public/clemenscords/multi_kmer/raw.csv --benchmark_out_format=csv --benchmark_repetitions=3 --benchmark_report_aggregates_only=false
int main(int argc, char** argv)
{

    register_all();

    benchmark::Initialize(&argc, argv);
    benchmark::RunSpecifiedBenchmarks();

    cleanup_csv("/srv/public/clemenscords/multi_kmer/final_raw.csv");
}