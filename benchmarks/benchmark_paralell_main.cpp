//
// Created by clem on 4/29/20.
// Copyright (c) 2020 Clemens Cords. All rights reserved.
//

#include <benchmarks/benchmark.hpp>
#include <seqan3/alphabet/nucleotide/dna4.hpp>

int main(int argc, char** argv)
{

    auto config = benchmark_arguments(100, 6000, 10000);

    benchmark::RegisterBenchmark("multi_par", &multi_kmer_construction<seqan3::dna4, true, 3, 4, 5, 6, 7, 8, 9, 10>, config, 8);
    benchmark::RegisterBenchmark("multi_seq", &multi_kmer_construction<seqan3::dna4, true, 3, 4, 5, 6, 7, 8, 9, 10>, config, 1);

    benchmark::Initialize(&argc, argv);
    benchmark::RunSpecifiedBenchmarks();
}