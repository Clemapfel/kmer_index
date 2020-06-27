
//
// Created by clem on 6/9/20.
// Copyright (c) 2020 Clemens Cords. All rights reserved.
//

#include <type_traits>
#include <cstdint>
#include <cassert>
#include <cmath>
#include <algorithm>
#include <unordered_map>

#include <fast_pow.hpp>
#include <benchmarks/input_generator.hpp>
#include <seqan3/alphabet/nucleotide/dna4.hpp>
#include <benchmark/benchmark.h>
#include <benchmarks/cleanup_csv.hpp>

using alphabet_t = seqan3::dna4;

template<size_t k>
struct kmer_hash
{
    private:
        inline static size_t _sigma = seqan3::alphabet_size<alphabet_t>;

        template<typename iterator_t>
        static size_t hash_aux_aux(iterator_t query_it, size_t i)
        {
            return seqan3::to_rank(*query_it) * kmer::detail::fast_pow(_sigma, k - i - 1);
        }

        template<typename iterator_t, size_t... is>
        static size_t hash_aux(iterator_t query_it, std::index_sequence<is...> sequence)
        {
            return (... + hash_aux_aux(query_it++, is));
        }

    public:
        // regular hash
        static size_t hash_for(std::vector<alphabet_t>& query)
        {
            auto it = query.begin();
            size_t size = k;

            size_t hash = 0;
            for (size_t i = 0; i < size; ++i)
                hash += (seqan3::to_rank(*it++) * kmer::detail::fast_pow(_sigma, size - i - 1));

            return hash;
        }

        // fold hash
        static size_t hash_fold(std::vector<alphabet_t>& query)
        {
            auto it = query.begin();
            return hash_aux(it, std::make_index_sequence<k>());
        }
};

size_t seed = 1234;
const size_t n_queries = 1000000;

// hash done with for
template<size_t k>
static void hash_for(benchmark::State& state)
{
    state.counters["seed"] = seed;

    input_generator<alphabet_t> input(seed++);
    std::vector<std::vector<alphabet_t>> kmers;
    kmers.reserve(n_queries);

    for (size_t i = 0; i < n_queries; ++i)
        kmers.push_back(input.generate_sequence(k));

    size_t i = 0;
    for (auto _ : state)
    {
        benchmark::DoNotOptimize(kmer_hash<k>::hash_for(kmers.at(i)));
        i = (i + 1 < kmers.size() ? i+1 : 0);
    }

    state.counters["k"] = k;
}

// hash done with for
template<size_t k>
static void hash_fold(benchmark::State& state)
{
    state.counters["seed"] = seed;

    input_generator<alphabet_t> input(seed++);
    std::vector<std::vector<alphabet_t>> kmers;
    kmers.reserve(n_queries);

    for (size_t i = 0; i < n_queries; ++i)
        kmers.push_back(input.generate_sequence(k));

    size_t i = 0;
    for (auto _ : state)
    {
        benchmark::DoNotOptimize(kmer_hash<k>::hash_fold(kmers.at(i)));
        i = (i + 1 < kmers.size() ? i+1 : 0);
    }

    state.counters["k"] = k;
}

template<size_t... ks>
void register_all()
{
    (benchmark::RegisterBenchmark("hash_fold", hash_fold<ks>), ...);
    (benchmark::RegisterBenchmark("hash_for", hash_for<ks>), ...);
    (benchmark::RegisterBenchmark("hash_seqan", hash_seqan<ks>), ...);
}

int main(int argc, char** argv)
{
    register_all<5, 10, 15, 20, 25>();

    benchmark::Initialize(&argc, argv);
    benchmark::RunSpecifiedBenchmarks();

    cleanup_csv("/home/clem/Documents/Workspace/kmer_index/source/benchmarks/hash_vs_hash/raw.csv");
}


