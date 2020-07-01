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
#include <seqan3/search/search.hpp>

#include <memory>

// benchmark intended to show that for looking up kmers, kmer-index is better than fm

size_t n_benchmarks_registered = 0;
size_t seed = 200;

using alphabet_t = seqan3::dna4;


// [SEQ] exact fm search
template<typename index_t>
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

// [SEQ] kmer exact search
template<size_t k>
static void kmer_search(benchmark::State& state, size_t query_length, const std::vector<alphabet_t>& text)
{
    state.counters["seed"] = seed;

    input_generator<alphabet_t> input(seed);
    auto index = kmer::make_kmer_index<k>(text);

    size_t i = 0;
    for (auto _ : state)
    {
        benchmark::DoNotOptimize(index.search(input.generate_sequence(query_length)));
    }

    state.counters["text_length"] = text.size();
    state.counters["query_length"] = query_length;
    state.counters["alphabet_size"] = seqan3::alphabet_size<alphabet_t>;

    seed++;
}

template<size_t... ks, typename index_t>
void register_all(std::vector<alphabet_t>& text, index_t* fm)
{
    (benchmark::RegisterBenchmark("kmer", &kmer_search<ks>, ks, text), ...);
    (benchmark::RegisterBenchmark("fm", &fm_search<index_t>, ks, text.size(), fm), ...);
}

//nohup ./JUST_K_BENCHMARK --benchmark_format=console --benchmark_counters_tabular=true --benchmark_out=/srv/public/clemenscords/just_k/1e7_to_1e10_raw.csv --benchmark_out_format=csv --benchmark_repetitions=15 --benchmark_report_aggregates_only=false

int main(int argc, char** argv)
{
    // 1
    auto input = input_generator<alphabet_t>(seed);
    auto text_1 = input.generate_sequence(1e8);
    auto fm_1 = seqan3::fm_index(text_1);

    register_all<3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24, 25, 26, 27, 28, 29, 30>(
            text_1, &fm_1);

    // 2
    auto text_2 = input.generate_sequence(1e7);
    auto fm_2 = seqan3::fm_index(text_2);

    register_all<3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24, 25, 26, 27, 28, 29, 30>(
            text_2, &fm_2);

    benchmark::Initialize(&argc, argv);
    benchmark::RunSpecifiedBenchmarks();
}