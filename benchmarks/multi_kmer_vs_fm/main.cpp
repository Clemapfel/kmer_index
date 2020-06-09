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
#include <seqan3/search/fm_index/bi_fm_index.hpp>
#include <seqan3/search/search.hpp>

size_t n_benchmarks_registered = 0;
size_t seed = 200;

// [SEQ] exact fm search
template<seqan3::alphabet alphabet_t>
static void fm_search(benchmark::State& state, size_t text_length, size_t query_length)
{
    state.counters["seed"] = seed;

    bool error_occurred = false;
    try
    {
        input_generator<alphabet_t> input(seed);
        auto text = input.generate_sequence(text_length);
        auto queries = input.generate_queries(100000, query_length);

        auto index = seqan3::fm_index(text);

        size_t i = 0;
        for (auto _ : state)
        {
            benchmark::DoNotOptimize(seqan3::search(queries.at(i), index));
            i = i < queries.size() - 1 ? i + 1 : 0;
        }
    }
    catch (...)
    {
        error_occurred = true;
    }

    state.counters["text_length"] = text_length;
    state.counters["query_length"] = query_length;
    state.counters["alphabet_size"] = seqan3::alphabet_size<alphabet_t>;
    state.counters["error_occurred"] = error_occurred;

    seed++;
}

// [SEQ] multi kmer search
template<seqan3::alphabet alphabet_t, size_t... ks>
static void multi_kmer_search(benchmark::State& state, size_t text_length, size_t query_length)
{
    state.counters["seed"] = seed;

    bool error_occurred = false;
    try
    {
        input_generator<alphabet_t> input(seed);
        auto text = input.generate_sequence(text_length);
        auto queries = input.generate_queries(100000, query_length);

        auto index = kmer::make_kmer_index<ks...>(text);

        size_t i = 0;
        for (auto _ : state)
        {
            benchmark::DoNotOptimize(index.search(queries.at(i)));
            i = i < queries.size() - 1 ? i + 1 : 0;
        }
    }
    catch (...)
    {
        error_occurred = true;
    }

    state.counters["text_length"] = text_length;
    state.counters["query_length"] = query_length;
    state.counters["alphabet_size"] = seqan3::alphabet_size<alphabet_t>;
    state.counters["error_occurred"] = error_occurred;

    std::vector<size_t> _all_ks = {ks...};
    size_t i = 0;
    for (auto k : _all_ks)
    {
        state.counters["k_" + std::to_string(i)] = k;
        i++;
    }

    seed++;
}

constexpr size_t text_length = 2000000;

template<size_t... ks>
void register_all()
{
    for (size_t query_length = 130; query_length <= 170; ++query_length)
    {

        benchmark::RegisterBenchmark("multi_kmer_search", &multi_kmer_search<seqan3::dna4, ks...>, text_length, query_length);
        benchmark::RegisterBenchmark("fm_search", &fm_search<seqan3::dna4>, text_length, query_length);

        n_benchmarks_registered += 2;
    }

    seqan3::debug_stream << n_benchmarks_registered << " benchmarks registered.\n";
}

//./KMER_VS_FM_BENCHMARK --benchmark_format=console --benchmark_counters_tabular=true --benchmark_out=/srv/public/clemenscords/kmer_vs_fm/raw.csv --benchmark_out_format=csv --benchmark_repetitions=10 --benchmark_report_aggregates_only=false

int main(int argc, char** argv)
{
    register_all<7, 10, 13, 17, 21, 23>();

    benchmark::Initialize(&argc, argv);
    benchmark::RunSpecifiedBenchmarks();

    //cleanup_csv("/srv/public/clemenscords/kmer_vs_fm/raw.csv");
}
