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

//nohup ./JUST_K_BENCHMARK --benchmark_format=console --benchmark_counters_tabular=true --benchmark_out=/srv/public/clemenscords/just_k/to_1e7_complere_raw.csv --benchmark_out_format=csv --benchmark_repetitions=100 --benchmark_report_aggregates_only=false

int main(int argc, char** argv)
{
    size_t text_length = 1e8;

    auto input = input_generator<seqan3::dna4>(seed);
    auto text = input.generate_sequence(text_length);
    auto fm = seqan3::fm_index{text};
    auto bi_fm = seqan3::bi_fm_index{text};

    for (size_t query_length = 3; query_length < 300; ++query_length)
    {
        benchmark::RegisterBenchmark("fm", &fm_search<decltype(fm)>, query_length, text_length, &fm);
        benchmark::RegisterBenchmark("bi_fm", &fm_search<decltype(bi_fm)>, query_length, text_length, &bi_fm);
    }
    benchmark::Initialize(&argc, argv);
    benchmark::RunSpecifiedBenchmarks();
}