//
// Created by clem on 6/20/20.
// Copyright (c) 2020 Clemens Cords. All rights reserved.
//
#include <kmer_index.hpp>
#include <benchmarks/input_generator.hpp>
#include <fast_pow.hpp>
#include <benchmarks/cleanup_csv.hpp>

#include "../benchmark/include/benchmark/benchmark.h"
#include <seqan3/alphabet/nucleotide/dna4.hpp>
#include <seqan3/search/fm_index/fm_index.hpp>
#include <seqan3/search/search.hpp>

#include <memory>
#include <seqan3/alphabet/nucleotide/dna5.hpp>
#include <seqan3/alphabet/nucleotide/dna15.hpp>

// benchmark intended to show that for looking up kmers, kmer-index is better than fm

size_t n_benchmarks_registered = 0;
size_t seed = 200;

// [SEQ] kmer exact search
template<seqan3::alphabet alphabet_t, size_t k>
static void kmer_search(benchmark::State& state, const std::vector<alphabet_t>& text)
{
    state.counters["seed"] = seed;

    input_generator<alphabet_t> input(seed);
    auto index = kmer::make_kmer_index<k>(text);

    size_t i = 0;
    for (auto _ : state)
    {
        auto query = input.generate_sequence(k);
        benchmark::DoNotOptimize(index.search(query));
    }

    state.counters["text_length"] = text.size();
    state.counters["query_length"] = k;
    state.counters["alphabet_size"] = seqan3::alphabet_size<alphabet_t>;

    seed++;
}

template<typename alphabet_t, size_t... ks>
void register_all(std::vector<alphabet_t>& text)
{
    (benchmark::RegisterBenchmark("kmer", &kmer_search<alphabet_t, ks>, text), ...);
}

constexpr size_t text_length = 1e6;

int main(int argc, char** argv)
{
    using namespace seqan3;

    // dna4
    auto input_4 = input_generator<dna4>(seed++);
    auto dna4_text = input_4.generate_sequence(text_length);
    register_all<dna4, 5, 10, 15>(dna4_text);

    // dna5
    auto input_5 = input_generator<dna5>(seed++);
    auto dna5_text = input_5.generate_sequence(text_length);
    register_all<dna5, 5, 10, 15>(dna5_text);

    // dna15
    auto input_15 = input_generator<dna15>(seed++);
    auto dna15_text = input_15.generate_sequence(text_length);
    register_all<dna15, 5, 10, 15>(dna15_text);

    benchmark::Initialize(&argc, argv);
    benchmark::RunSpecifiedBenchmarks();
}