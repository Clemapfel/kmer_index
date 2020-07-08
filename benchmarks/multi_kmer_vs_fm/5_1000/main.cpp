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

// [SEQ] multi kmer exact search
template<typename index_t>
static void kmer_search(benchmark::State& state, size_t query_length, size_t text_length, const index_t* index)
{
    state.counters["seed"] = seed;

    input_generator<alphabet_t> input(seed);

    size_t i = 0;
    for (auto _ : state)
    {
        auto query = input.generate_sequence(query_length);
        benchmark::DoNotOptimize(index->search(query));
    }

    state.counters["text_length"] = text_length;
    state.counters["query_length"] = query_length;
    state.counters["alphabet_size"] = seqan3::alphabet_size<alphabet_t>;

    seed++;
}

template<typename fm_index_t, typename kmer_index_t>
void register_all(size_t query_length, size_t text_length, fm_index_t* fm, kmer_index_t* kmer)
{
    benchmark::RegisterBenchmark("kmer", &kmer_search<kmer_index_t>, query_length, text_length, kmer);
    benchmark::RegisterBenchmark("fm", &fm_search<fm_index_t>, query_length, text_length, fm);
}

constexpr size_t text_length = 1e4;

int main(int argc, char** argv)
{
    auto input = input_generator<alphabet_t>(seed);
    auto text = input.generate_sequence(text_length);

    auto fm = seqan3::fm_index(text);
    auto kmer = kmer::make_kmer_index<4, 5, 6, 7, 8, 9, 10, 11, 13, 15, 17, 19, 21, 23, 25, 27, 29, 31>(text, std::thread::hardware_concurrency());

    for (size_t i = 32; i <= 1000; ++i)
        register_all<decltype(fm), decltype(kmer)>(i, text_length, &fm, &kmer);

    benchmark::Initialize(&argc, argv);
    benchmark::RunSpecifiedBenchmarks();
}