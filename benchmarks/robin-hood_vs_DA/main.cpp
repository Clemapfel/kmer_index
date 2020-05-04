//
// Created by clem on 5/4/20.
// Copyright (c) 2020 Clemens Cords. All rights reserved.
//

#include "../benchmark.hpp"

// test runtime for DA vs robin-hood map

template<seqan3::alphabet alphabet_t, size_t k>
void register_benchmarks()
{
    // setup configs
    std::vector<benchmark_arguments<alphabet_t>> configs;

    for (size_t text_size : {10000, 100000, 1000000})
    {
        for (size_t query_size : {k, 10 * k})
        {
            configs.push_back(benchmark_arguments<alphabet_t>(query_size, 100000, text_size));
        }
    }

    // register benchmarks
    for (auto config : configs)
    {
        benchmark::RegisterBenchmark("robin_Hood", &kmer_construction<alphabet_t, false, k>, config);
        benchmark::RegisterBenchmark("DA", &kmer_construction<alphabet_t, true, k>, config);
    }
}

template<size_t... ks>
void register_all(int argc, char** argv)
{
    (register_benchmarks<seqan3::dna4, ks>(), ...);

    benchmark::Initialize(&argc, argv);
    benchmark::RunSpecifiedBenchmarks();
}

int main(int argc, char** argv)
{
    register_all<5, 6, 7, 8, 9, 10>(argc, argv);
}