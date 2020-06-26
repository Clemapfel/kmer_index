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

// fm
template<seqan3::alphabet alphabet_t>
static void fm_search(benchmark::State& state, size_t query_length, size_t text_length)
{
    state.counters["seed"] = seed;

    auto input = input_generator<seqan3::dna4>(seed);
    auto text = input.generate_sequence(text_length);
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

// kmer
template<seqan3::alphabet alphabet_t, size_t k>
static void kmer_search(benchmark::State& state, size_t query_length, size_t text_length)
{
    state.counters["seed"] = seed;
    auto input = input_generator<seqan3::dna4>(seed);
    auto text = input.generate_sequence(text_length);
    auto index = kmer::make_kmer_index<k>(text);

    bool error_occurred = false;
    try
    {
        size_t i = 0;
        for (auto _ : state)
        {
            benchmark::DoNotOptimize(index.search(input.generate_sequence(query_length)));
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
    state.counters["valid"] = not error_occurred;

    seed++;
}

template<size_t k>
void register_kmer(size_t text_length)
{
    benchmark::RegisterBenchmark("kmer_single", &kmer_search<seqan3::dna4, k>, k, text_length);
}


template<size_t... ks>
void register_all()
{
    for (size_t text_length : {1e9, 1e10})
    {
        (register_kmer<ks>(text_length), ...);
        (benchmark::RegisterBenchmark("single_kmer", &fm_search<seqan3::dna4>, ks, text_length), ...);
    }
}

int main(int argc, char** argv)
{
    register_all<3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24, 25, 26, 27, 28, 29, 30>();

    benchmark::Initialize(&argc, argv);
    benchmark::RunSpecifiedBenchmarks();
}