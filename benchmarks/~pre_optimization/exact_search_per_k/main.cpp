//
// Created by clem on 5/4/20.
// Copyright (c) 2020 Clemens Cords. All rights reserved.
//

#include <benchmarks/cleanup_csv.hpp>
#include "../benchmark.hpp"

inline std::atomic<size_t> N_BENCHMARKS = 0;

template<seqan3::alphabet alphabet_t, size_t k>
void register_benchmarks()
{
    size_t text_size = 500000;

    int int_k = int(k);
    for (int offset = -10; offset <= 3*int_k; offset++)
    {
        if (int_k + offset < k-2 or int_k > 200)
            continue;

        auto config = benchmark_arguments<alphabet_t>(int_k + offset, 100000, text_size);

        benchmark::RegisterBenchmark("kmer_search_map", &kmer_search<alphabet_t, false, k>, config);
        if (k <= 10)
            benchmark::RegisterBenchmark("kmer_search_da", &kmer_search<alphabet_t, true, k>, config);
        benchmark::RegisterBenchmark("fm_search", &fm_search<alphabet_t>, config);
        N_BENCHMARKS++;
    }

    seqan3::debug_stream << std::to_string(N_BENCHMARKS) + " benchmarks currently registered\n";
}

template<size_t... ks>
void register_all(int argc, char** argv)
{
    (register_benchmarks<seqan3::dna4, ks>(), ...);

    benchmark::Initialize(&argc, argv);
    benchmark::RunSpecifiedBenchmarks();
}

/* EXECUTABLE ARGUMENTS:
 * --benchmark_out=/home/clem/Documents/Workspace/kmer_index/source/benchmarks/exact_search_per_k/raw.csv

 *
--benchmark_format=console
--benchmark_counters_tabular=true
--benchmark_out=/srv/public/clemenscords/exact_search_per_k/raw.csv
--benchmark_out_format=csv
--benchmark_repetitions=10
--benchmark_report_aggregates_only=true
 */

int main(int argc, char** argv)
{
    register_all<13, 14, 15, 16, 17, 18, 19, 20>(argc, argv);
    //cleanup_csv("/home/clem/Documents/Workspace/kmer_index/source/benchmarks/exact_search_per_k/raw.csv");
}