//
// Created by clem on 5/4/20.
// Copyright (c) 2020 Clemens Cords. All rights reserved.
//

#include <benchmarks/cleanup_csv.hpp>
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
        //benchmark::RegisterBenchmark("robin_hood_construction", &kmer_construction<alphabet_t, false, k>, config);
        //benchmark::RegisterBenchmark("DA_construction", &kmer_construction<alphabet_t, true, k>, config);
        benchmark::RegisterBenchmark("robin_hood_search", &kmer_search<alphabet_t, false, k>, config);
        benchmark::RegisterBenchmark("DA_search", &kmer_search<alphabet_t, true, k>, config);
    }
}

template<size_t... ks>
void register_all(int argc, char** argv)
{
    (register_benchmarks<seqan3::dna4, ks>(), ...);

    benchmark::Initialize(&argc, argv);
    benchmark::RunSpecifiedBenchmarks();
}

void set_local_env()
{
    std::string benchmark_out = "BENCHMARK_OUT=/home/clem/Documents/Workspace/kmer_index/source/benchmarks/robin-hood_vs_DA/test.csv";
    std::string out_format = "BENCHMARK_OUT_FORMAT=csv";
    std::string console_format = "BENCHMARK_FORMAT=csv";
    std::string tabular = "BENCHMARK_COUNTERS_TABULAR=true";

    putenv(benchmark_out.data());
    putenv(out_format.data());
    putenv(console_format.data());
    putenv(tabular.data());

}

/* EXECUTABLE ARGUMENTS:
--benchmark_format=console
--benchmark_counters_tabular=true
--benchmark_out=/home/clem/Documents/Workspace/kmer_index/source/benchmarks/robin-hood_vs_DA/raw.csv"
--benchmark_out=
--benchmark_out_format=csv
--benchmark_min_time=0.5
--benchmark_repetitions=1
--benchmark_report_aggregates_only=false
 */

int main(int argc, char** argv)
{
    register_all<5, 6, 7, 8, 9, 10, 15, 20>(argc, argv);
    //cleanup_csv("/home/clem/Documents/Workspace/kmer_index/source/benchmarks/robin-hood_vs_DA/raw.csv");
}