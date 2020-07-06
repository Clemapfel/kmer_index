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
        auto query = input.generate_sequence(query_length);
        benchmark::DoNotOptimize(search(query, *index));
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
        auto query = input.generate_sequence(query_length);
        benchmark::DoNotOptimize(index.search(query));
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

template<size_t... ks>
void register_kmer(std::vector<alphabet_t>& text)
{
    (benchmark::RegisterBenchmark("kmer", &kmer_search<ks>, ks, text), ...);
}

//nohup ./JUST_K_BENCHMARK --benchmark_format=console --benchmark_counters_tabular=true --benchmark_out=/srv/public/clemenscords/just_k/to_1e7_complere_raw.csv --benchmark_out_format=csv --benchmark_repetitions=100 --benchmark_report_aggregates_only=false

int main(int argc, char** argv)
{
    auto input = input_generator<alphabet_t>(seed);

    /*
    auto text_1 = input.generate_sequence(1e3);
    auto fm_1 = seqan3::fm_index(text_1);

    register_all<3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24, 25, 26, 27, 28, 29, 30>(
            text_1, &fm_1);

    auto text_2 = input.generate_sequence(1e4);
    auto fm_2 = seqan3::fm_index(text_2);

    register_all<3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24, 25, 26, 27, 28, 29, 30>(
            text_2, &fm_2);

    auto text_3 = input.generate_sequence(1e5);
    auto fm_3 = seqan3::fm_index(text_3);

    register_all<3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24, 25, 26, 27, 28, 29, 30>(
            text_3, &fm_3);

    auto text_4 = input.generate_sequence(1e6);
    auto fm_4 = seqan3::fm_index(text_4);

    register_all<3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24, 25, 26, 27, 28, 29, 30>(
            text_4, &fm_4);

    auto text_5 = input.generate_sequence(1e7);
    auto fm_5 = seqan3::fm_index(text_5);

    register_all<3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24, 25, 26, 27, 28, 29, 30>(
            text_5, &fm_5);

    auto text_6 = input.generate_sequence(1e8);
    auto fm_6 = seqan3::fm_index(text_6);

    register_all<3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24, 25, 26, 27, 28, 29, 30>(
            text_6, &fm_6);
    */
    /*
    auto text_7 = input.generate_sequence(1e9);
    auto fm_7 = seqan3::fm_index(text_7);

    register_all<3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24, 25, 26, 27, 28, 29, 30>(
            text_7, &fm_7);
            */

    auto text_8 = input.generate_sequence(1e8);
    register_kmer<3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24, 25, 26, 27, 28, 29, 30>(
            text_8);

    benchmark::Initialize(&argc, argv);
    benchmark::RunSpecifiedBenchmarks();
}