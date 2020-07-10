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

size_t seed = 200;

// [SEQ] seqan3 fm search
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

// [SEQ] kmer search
template<seqan3::alphabet alphabet_t, typename index_t>
static void multi_kmer_search(benchmark::State& state, size_t query_length, size_t text_length, index_t* index)
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

    std::vector<size_t> all_ks = index->get_ks();
    size_t j = 0;
    for (auto k : all_ks)
    {
        state.counters["k_" + std::to_string(j)] = k;
        j++;
    }
    state.counters["valid"] = not error_occurred;

    seed++;
}

// benchmark shows that one multi_kmer can handle more than just the exact k given
Kar= 1e6;

int main(int argc, char** argv)
{
    auto input = input_generator<seqan3::dna4>(seed);
    auto text = input.generate_sequence(text_length);

    auto multi_kmer = kmer::make_kmer_index<5, 7, 9, 11, 13, 15, 17, 19, 21, 23, 25, 27, 29>(text);
    auto single_kmer = kmer::make_kmer_index<10>(text);
    auto fm = seqan3::fm_index(text);

    for (size_t query_length = 4; query_length <= 50; ++query_length)
    {
        benchmark::RegisterBenchmark("multi_kmer", &multi_kmer_search<seqan3::dna4, decltype(multi_kmer)>, query_length, text_length, &multi_kmer);
        benchmark::RegisterBenchmark("single_kmer", &multi_kmer_search<seqan3::dna4, decltype(single_kmer)>, query_length, text_length, &single_kmer);
        benchmark::RegisterBenchmark("fm", &fm_search<seqan3::dna4, decltype(fm)>, query_length, text_length, &fm);
    }

    benchmark::Initialize(&argc, argv);
    benchmark::RunSpecifiedBenchmarks();
}