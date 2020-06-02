//
// Created by clem on 6/1/20.
// Copyright (c) 2020 Clemens Cords. All rights reserved.
//

//
// Created by clem on 5/30/20.
// Copyright (c) 2020 Clemens Cords. All rights reserved.
//
#include <kmer_index.hpp>
#include <benchmarks/input_generator.hpp>
#include <fast_pow.hpp>
#include <benchmarks/cleanup_csv.hpp>

#include <benchmark/benchmark.h>

#include <seqan3/alphabet/nucleotide/dna4.hpp>
#include <seqan3/search/fm_index/fm_index.hpp>
#include <seqan3/search/search.hpp>

size_t n_benchmarks_registered = 0;
size_t seed = 0;

// [SEQ] fm construction
template<seqan3::alphabet alphabet_t>
static void fm_construction(benchmark::State& state, size_t text_length)
{
    input_generator<alphabet_t> input(seed++);
    auto text = input.generate_sequence(text_length);

    auto index = seqan3::fm_index(text);

    size_t i = 0;
    for (auto _ : state)
    {
        benchmark::DoNotOptimize(seqan3::fm_index(text));
    }

    state.counters["text_length"] = text_length;
    state.counters["alphabet_size"] = seqan3::alphabet_size<alphabet_t>;
    state.counters["n_threads"] = 1;
}

// [SEQ] single kmer construction
template<seqan3::alphabet alphabet_t, size_t k>
static void kmer_construction(benchmark::State& state, size_t text_length)
{
    input_generator<alphabet_t> input(seed++);
    auto text = input.generate_sequence(text_length);

    size_t i = 0;
    for (auto _ : state)
    {
        benchmark::DoNotOptimize(kmer::make_kmer_index<k>(text));
    }

    state.counters["text_length"] = text_length;
    state.counters["alphabet_size"] = seqan3::alphabet_size<alphabet_t>;
    state.counters["k_0"] = k;
    state.counters["n_threads"] = 1;
}

// [SEQ] multi kmer construction
template<seqan3::alphabet alphabet_t, size_t... ks>
static void seq_multi_kmer_construction(benchmark::State& state, size_t text_length)
{
    input_generator<alphabet_t> input(seed++);
    auto text = input.generate_sequence(text_length);

    size_t i = 0;
    for (auto _ : state)
    {
        benchmark::DoNotOptimize(kmer::make_kmer_index<ks...>(text, 1));
    }

    state.counters["text_length"] = text_length;
    state.counters["alphabet_size"] = seqan3::alphabet_size<alphabet_t>;

    size_t k_i = 0;
    for (size_t k : std::vector<size_t>{(ks, ...)})
    {
        state.counters["k_" + std::to_string(k_i++)] = k;
    }
    state.counters["n_threads"] = 1;
}

// [PAR] multi kmer construction
template<seqan3::alphabet alphabet_t, size_t... ks>
static void par_multi_kmer_construction(benchmark::State& state, size_t text_length)
{
    input_generator<alphabet_t> input(seed++);
    auto text = input.generate_sequence(text_length);

    std::vector<size_t> all_ks{(ks, ...)};
    size_t n_threads = std::min<size_t>(all_ks.size(), std::thread::hardware_concurrency());

    size_t i = 0;
    for (auto _ : state)
    {
        benchmark::DoNotOptimize(kmer::make_kmer_index<ks...>(text, n_threads));
    }

    state.counters["text_length"] = text_length;
    state.counters["alphabet_size"] = seqan3::alphabet_size<alphabet_t>;

    size_t k_i = 0;
    for (size_t k : all_ks)
    {
        state.counters["k_" + std::to_string(k_i++)] = k;
    }
    state.counters["n_threads"] = n_threads;
}

template<size_t k>
void register_all(size_t min_text_length, size_t max_text_length, size_t step_size)
{
    constexpr size_t k_0 = k - 1;
    constexpr size_t k_1 = k;
    constexpr size_t k_2 = k + 1;

    for (size_t text_length = min_text_length; text_length <= max_text_length; text_length += step_size)
    {
        benchmark::RegisterBenchmark()
    }
}

//./KMER_VS_FM_BENCHMARK --benchmark_format=console --benchmark_counters_tabular=true --benchmark_out=/srv/public/clemenscords/main_benchmark/raw.csv --benchmark_out_format=csv --benchmark_repetitions=100 --benchmark_report_aggregates_only=false

int main(int argc, char** argv)
{
    TODO: what to show?

    register_all<1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20>();

    benchmark::Initialize(&argc, argv);
    benchmark::RunSpecifiedBenchmarks();

    cleanup_csv("/srv/public/clemenscords/main_benchmark/raw.csv");
}


