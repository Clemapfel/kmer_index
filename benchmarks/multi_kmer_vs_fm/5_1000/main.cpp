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

    for (auto _ : state)
    {
        auto query = input.generate_sequence(query_length);
        benchmark::DoNotOptimize(index->search(query));
    }

    state.counters["text_length"] = text_length;
    state.counters["query_length"] = query_length;
    state.counters["alphabet_size"] = seqan3::alphabet_size<alphabet_t>;

    size_t i = 0;
    for (auto k : index->get_ks())
    {
        std::string label = "k_" + std::to_string(i);
        state.counters[label] = k;
    }

    seed++;
}

// [SEQ] multi kmer exact search
template<typename index_t>
static void kmer_experimental_search(benchmark::State& state, size_t query_length, size_t text_length, const index_t* index)
{
    state.counters["seed"] = seed;

    input_generator<alphabet_t> input(seed);

    for (auto _ : state)
    {
        auto query = input.generate_sequence(query_length);
        benchmark::DoNotOptimize(index->experimental_search(query));
    }

    state.counters["text_length"] = text_length;
    state.counters["query_length"] = query_length;
    state.counters["alphabet_size"] = seqan3::alphabet_size<alphabet_t>;

    size_t i = 0;
    for (auto k : index->get_ks())
    {
        std::string label = "k_" + std::to_string(i);
        state.counters[label] = k;
    }

    seed++;
}

template<typename fm_index_t, typename kmer_index_t>
void register_all(size_t query_length, size_t text_length, fm_index_t* fm, kmer_index_t* kmer)
{
    benchmark::RegisterBenchmark("kmer", &kmer_search<kmer_index_t>, query_length, text_length, kmer);
    benchmark::RegisterBenchmark("kmer_experimental", &kmer_experimental_search<kmer_index_t>, query_length, text_length, kmer);
    benchmark::RegisterBenchmark("fm", &fm_search<fm_index_t>, query_length, text_length, fm);
}

constexpr size_t text_length = 1e5;

int main(int argc, char** argv)
{
    auto input = input_generator<alphabet_t>(seed);
    auto text = input.generate_sequence(text_length);

    auto fm = seqan3::fm_index(text);
    auto kmer = kmer::make_kmer_index<5, 6, 7, 8, 9, 10, 11, 12, 13, 15, 17, 19, 21, 23, 25, 27, 29, 31>(text, std::thread::hardware_concurrency());

    std::vector<size_t> all_primes = {19, 23, 29, 31, 37, 41, 43, 47, 53, 59, 61, 67,
 71, 73, 79, 83, 89, 97, 101, 103, 107, 109, 113, 127, 131, 137, 139, 149,
151, 157, 163, 167, 173, 179, 181, 191, 193, 197, 199, 211, 223, 227, 229,
 233, 239, 241, 251, 257, 263, 269, 271, 277, 281, 283, 293, 307, 311,
313, 317, 331, 337, 347, 349, 353, 359, 367, 373, 379, 383, 389, 397, 401,
 409, 419, 421, 431, 433, 439, 443, 449, 457, 461, 463, 467, 479, 487,
491, 499, 503, 509, 521, 523, 541, 547, 557, 563, 569, 571, 577, 587, 593,
599, 601, 607, 613, 617, 619, 631, 641, 643, 647, 653, 659, 661, 673, 677,
 683, 691, 701, 709, 719, 727, 733, 739, 743, 751, 757, 761, 769, 773,
787, 797, 809, 811, 821, 823, 827, 829, 839, 853, 857, 859, 863, 877, 881,
 883, 887, 907, 911, 919, 929, 937, 941, 947, 953, 967, 971, 977, 983, 991, 997};

    for (size_t p : all_primes)
    {
        register_all<decltype(fm), decltype(kmer)>(p-1, text_length, &fm, &kmer);
        register_all<decltype(fm), decltype(kmer)>(p, text_length, &fm, &kmer);
        register_all<decltype(fm), decltype(kmer)>(p+1, text_length, &fm, &kmer);
    }

    benchmark::Initialize(&argc, argv);
    benchmark::RunSpecifiedBenchmarks();
}