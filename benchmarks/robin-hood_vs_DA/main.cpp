//
// Created by clem on 5/4/20.
// Copyright (c) 2020 Clemens Cords. All rights reserved.
//

#include "../benchmark.hpp"

// test runtime for DA vs robin-hood map

template<seqan3::alphabet alphabet_t, size_t... ks>
void register_all()
{
    
}

int main(int argc, char** argv)
{
    benchmark::RegisterBenchmark("robin_hood", &kmer_construction<seqan3::dna4, false, 5>, config);
    benchmark::RegisterBenchmark("da", &kmer_construction<seqan3::dna4, true, 5>, config);

    auto config = benchmark_arguments<seqan3::dna4>(100, 6000, 10000);

    benchmark::RegisterBenchmark("multi_par", &multi_kmer_construction<seqan3::dna4, true, 3, 4, 5, 6, 7, 8, 9, 10>, config, 8)->UseRealTime();
    benchmark::RegisterBenchmark("multi_seq", &multi_kmer_construction<seqan3::dna4, true, 3, 4, 5, 6, 7, 8, 9, 10>, config, 1)->UseRealTime();

    benchmark::Initialize(&argc, argv);
    benchmark::RunSpecifiedBenchmarks();
}