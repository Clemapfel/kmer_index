// Copyright (c) 2020 Clemens Cords. All rights reserved.

#include <benchmarks/cleanup_csv.hpp>
#include <benchmarks/input_generator.hpp>
#include "../benchmark.hpp"

/* test if it's worth it to run heuristic before da or if reallocating DA table as it grows isn't that much slower */

// [SEQ] construction with da. with/without heuristic
template<seqan3::alphabet alphabet_t, bool use_heuristic, size_t k>
static void kmer_da_heuristic(benchmark::State& state, size_t text_length)
{
    debug::USE_DA_HEURISTIC = use_heuristic;

    auto gen = input_generator<alphabet_t>{1234};

    std::vector<std::vector<alphabet_t>> texts;
    for (size_t i = 0; i < 1000; ++i)
        texts.push_back(gen.generate_sequence(text_length));

    size_t i = 0;
    for (auto _ : state)
    {
        benchmark::DoNotOptimize(debug::make_kmer_index<k>(texts.at(i), true, 1));
        i = i+1 < texts.size() ? i+1 : 0;
    }

    state.counters["k"] = k;
    state.counters["text_size"] = text_length;
    state.counters["alphabet_size"] = seqan3::alphabet_size<alphabet_t>;
}

template<seqan3::alphabet alphabet_t, size_t k>
void register_benchmark()
{
    std::vector<benchmark_arguments<alphabet_t>> configs;

    for (size_t text_size : {10000, 100000, 500000})
    {
        benchmark::RegisterBenchmark("with_heuristic", &kmer_da_heuristic<alphabet_t, true, k>, text_size);
        benchmark::RegisterBenchmark("without_heuristic", &kmer_da_heuristic<alphabet_t, false, k>, text_size);
    }
}

template<size_t... ks>
void register_all(int argc, char** argv)
{
    (register_benchmark<seqan3::dna4, ks>(), ...);

    benchmark::Initialize(&argc, argv);
    benchmark::RunSpecifiedBenchmarks();
}

/* EXECUTABLE ARGUMENTS:
--benchmark_format=console
--benchmark_counters_tabular=true
--benchmark_out=/srv/public/clemenscords/main_benchmark/raw.csv
--benchmark_out_format=csv
--benchmark_repetitions=100
--benchmark_report_aggregates_only=false
 */

int main(int argc, char** argv)
{
    register_all<10, 11, 12, 13, 15>(argc, argv);
    cleanup_csv("/srv/public/clemenscords/main_benchmark/raw.csv");
}
