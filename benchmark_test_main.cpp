//
// Created by Clemens Cords on 3/3/20.
// Copyright (c) 2020 Clemens Cords. All rights reserved.
//

#include <benchmark.hpp>

/* EXECUTABLE ARGUMENTS:
--benchmark_format=console
--benchmark_counters_tabular=true
--benchmark_out=../source/benchmark_out_raw.csv
--benchmark_out_format=csv
--benchmark_min_time=0.5
--benchmark_repetitions=3
--benchmark_report_aggregates_only=false
 */

// main
int main(int argc, char** argv)
{
    std::vector<size_t> text_sizes = {1000, 10000);

    register_all_benchmarks<seqan3::dna4, false, 5, 6, 7>(text_sizes);
    register_all_benchmarks<seqan3::dna4, true, 5, 6, 7>(text_sizes);

    benchmark::Initialize(&argc, argv);
    benchmark::RunSpecifiedBenchmarks();

    cleanup_csv("../source/benchmark_out_raw.csv");
}

