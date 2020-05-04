//
// Created by Clemens Cords on 3/3/20.
// Copyright (c) 2020 Clemens Cords. All rights reserved.
//

#include <benchmarks/benchmark.hpp>

// benchmark: kmer_index construction based on
template<seqan3::alphabet alphabet_t, bool use_heuristic, size_t k>
static void kmer_construction_heuristic(benchmark::State& state, size_t text_size)
{
    auto text = input_generator<alphabet_t>::generate_text(text_size, {});

    detail::USE_DA_HEURISTIC = use_heuristic;

    for (auto _ : state)
        benchmark::DoNotOptimize(make_kmer_index<true, k>(text));

    state.counters["k"] = k;

    // log memory used
    auto index = make_kmer_index<true, k>(text);
    state.counters["memory_used(mb)"] = sizeof(index) / 1e6;
}

template<size_t... ks>
void register_all(size_t text_size)
{
    ((benchmark::RegisterBenchmark("heuristic", &kmer_construction_heuristic<seqan3::dna4, true, ks>, text_size),
     benchmark::RegisterBenchmark("no_heuristic", &kmer_construction_heuristic<seqan3::dna4, false, ks>, text_size)), ...);
}

// main
int main(int argc, char** argv)
{
    // --benchmark_format=console --benchmark_out=/srv/public/clemenscords/heuristic_out.csv --benchmark_out_format=csv --benchmark_repetitions=3 --benchmark_min_time=5

    size_t text_size = size_t(1e5);
    register_all<5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24, 25>(text_size);

    benchmark::Initialize(&argc, argv);
    benchmark::RunSpecifiedBenchmarks();

    std::cout << "finished running benchmarks.\n";
}

