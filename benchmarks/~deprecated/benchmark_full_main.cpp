//
// Created by Clemens Cords on 3/3/20.
// Copyright (c) 2020 Clemens Cords. All rights reserved.
//

#include "benchmark.hpp"

/* EXECUTABLE ARGUMENTS:
--benchmark_format=console
--benchmark_counters_tabular=true
--benchmark_out=/srv/public/clemenscords/benchmark_out_raw.csv
--benchmark_out_format=csv
--benchmark_min_time=0.5
--benchmark_repetitions=3
--benchmark_report_aggregates_only=false
 */

// main
int main(int argc, char** argv)
{
    std::vector<size_t> text_sizes = {1000, 10000, 100000, 1000000, 1000000000);

    register_all_benchmarks<seqan3::dna4, false, 5, 6, 7, 8, 9, 10, 15, 20>(text_sizes);
    register_all_benchmarks<seqan3::dna4, true, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24, 25, 26, 27, 28, 29, 30>(text_sizes);

    benchmark::Initialize(&argc, argv);
    benchmark::RunSpecifiedBenchmarks();

    std::cout << "finished running benchmarks.\n";
}

/*
// alphabets: dna4
// ks: 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24, 25, 26, 27, 28, 29, 30
// query_size_i = k_i -0.5*k, k_i, k+1 (-1 = average case, +1 = worst case);

// text_sizes: 3e9 = human dna, 249e6 = biggest chromosome,

// no_hit_ratio: 0, 1,

// native l1 cage = 4000 chars
// native l2 cage = 32000 chars
// native l3 cage = 768000

cache effect: 10k
    bacterielles genom: 10³, 10⁶, 10^9
    k: 8, log4(text_länge) = 10mer, genom = 16, 28

    dna4, 10er murphy

bitfields: speiceher mehrere postitoine von k in einem bit, häng A = 0 dran damit


*/




