//
// Created by clem on 4/29/20.
// Copyright (c) 2020 Clemens Cords. All rights reserved.
//

#include <benchmarks/benchmark.hpp>
#include <seqan3/alphabet/nucleotide/dna4.hpp>
#include <thread_pool.hpp>

static void test_benchmark(benchmark::State& state, std::string str)
{
    size_t i = 0;
    for (auto _ : state)
    {
        if (i < 3)
            //sync_print(str);
        ++i;
    }
}

int main(int argc, char** argv)
{
    auto config = benchmark_arguments<seqan3::dna4>(100, 6000, 10000);

    benchmark::RegisterBenchmark("multi_par", &multi_kmer_construction<seqan3::dna4, true, 3, 4, 5, 6, 7, 8, 9, 10>, config, 8)->UseRealTime();
    benchmark::RegisterBenchmark("multi_seq", &multi_kmer_construction<seqan3::dna4, true, 3, 4, 5, 6, 7, 8, 9, 10>, config, 1)->UseRealTime();

    benchmark::Initialize(&argc, argv);
    benchmark::RunSpecifiedBenchmarks();
}
