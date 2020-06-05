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
size_t seed = 0;

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
    state.counters["error_occured"] = error_occurred;

    seed++;
}


// [SEQ] bi fm search
template<seqan3::alphabet alphabet_t>
static void bi_fm_search(benchmark::State& state, size_t text_length, size_t query_length)
{
    state.counters["seed"] = seed;

    bool error_occurred = false;
    try
    {
        input_generator<alphabet_t> input(seed);
        auto text = input.generate_sequence(text_length);
        auto queries = input.generate_queries(100000, query_length);

        auto index = seqan3::bi_fm_index(text);

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
    state.counters["error_occured"] = error_occurred;

    seed++;
}


// [SEQ] exact single kmer search
template<seqan3::alphabet alphabet_t, size_t k>
static void kmer_search(benchmark::State& state, size_t text_length, size_t query_length)
{
    state.counters["seed"] = seed;

    bool error_occurred = false;
    try
    {
        input_generator<alphabet_t> input(seed);
        auto text = input.generate_sequence(text_length);
        auto queries = input.generate_queries(100000, query_length);

        auto index = kmer::make_kmer_index<k>(text);

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
    state.counters["k"] = k;
    state.counters["error_occured"] = error_occurred;
    seed++;
}

template<size_t k>
void register_kmer_benchmark(size_t text_length, size_t min_query_length, size_t max_query_length)
{
    for (size_t query_length = min_query_length; query_length < max_query_length; ++query_length)
    {
        benchmark::RegisterBenchmark("kmer_search", &kmer_search<seqan3::dna4, k>, text_length, query_length);
        n_benchmarks_registered += 1;
    }
}

constexpr size_t text_length = 1000000;

template<size_t... ks>
void register_all()
{
    (register_kmer_benchmark<ks>(text_length, ks - 2, 6*ks), ...);

    /*
    std::vector<size_t> _all_ks = {ks...};
    size_t min_k = std::numeric_limits<size_t>::max();
    size_t max_k = 0;

    for (size_t k : _all_ks)
    {
        if (k < min_k)
            min_k = k;
        if (k > max_k)
            max_k = k;
    }


    for (size_t query_length = 1; query_length < 6*max_k; query_length++)
    {
        if (query_length >= 149)
            benchmark::RegisterBenchmark("fm_search", &fm_search<seqan3::dna4>, text_length, query_length);

        benchmark::RegisterBenchmark("bi_fm_search", &fm_search<seqan3::dna4>, text_length, query_length);
        n_benchmarks_registered += 2;
    }*/

    seqan3::debug_stream << n_benchmarks_registered << " benchmarks registered.\n";
}

//./KMER_VS_FM_BENCHMARK --benchmark_format=console --benchmark_counters_tabular=true --benchmark_out=/srv/public/clemenscords/kmer_vs_fm/raw.csv --benchmark_out_format=csv --benchmark_repetitions=10 --benchmark_report_aggregates_only=false

int main(int argc, char** argv)
{
    register_all<17, 18, 19, 20, 21, 22, 23, 24, 25, 26, 27, 28, 29, 30>();

    benchmark::Initialize(&argc, argv);
    benchmark::RunSpecifiedBenchmarks();

    //cleanup_csv("/srv/public/clemenscords/kmer_vs_fm/raw.csv");

}
