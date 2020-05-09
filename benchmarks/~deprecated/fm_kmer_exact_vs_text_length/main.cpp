// Copyright (c) 2020 Clemens Cords. All rights reserved.

#include <benchmarks/input_generator.hpp>
#include <benchmarks/benchmark.hpp>
#include <benchmarks/cleanup_csv.hpp>

#include <random>
#include <vector>

// show kmer search time for steady query vs fm index

// generate queries used for simulation
template<seqan3::alphabet alphabet_t>
std::vector<std::vector<alphabet_t>> generate_query_set(size_t n_queries, float offlength_percentage = 0.05)
{
    size_t seed = 1234;

    size_t target_length = 150; // [1]
    auto gen = input_generator<alphabet_t>(seed);
    size_t n_offlength = n_queries * offlength_percentage; // number of reads that aren't 150

    std::vector<std::vector<alphabet_t>> output;

    // target lengths
    for (size_t i = 0; i < n_queries - n_offlength; ++i)
    {
        output.push_back(gen.generate_sequence(target_length));
    }

    // add in sequencing errors
    std::mt19937 engine(seed);
    std::uniform_int_distribution<int> dist(-2, 2);
    for (size_t i = 0; i < n_offlength; ++i)
    {
        output.push_back(gen.generate_sequence(target_length + dist(engine)));
    }

    return output;
}

// ###############################################################################

// [SEQ] KMER
template<seqan3::alphabet alphabet_t, size_t k>
static void kmer_search_vs_text_size(benchmark::State& state, size_t text_size)
{
    // generate text and queries
    auto gen = input_generator<alphabet_t>(1234);

    auto queries = generate_query_set<alphabet_t>(1e5);
    std::vector<alphabet_t> text = gen.generate_text(text_size, queries);

    auto index = debug::make_kmer_index<k, k+1, k+2>(text, false, 3); // equivalent to stddev

    size_t i = 0;
    for (auto _ : state)
    {
        benchmark::DoNotOptimize(index.search(queries.at(i)));
        i = i < queries.size()-1 ? i + 1 : 0;   // cycle through queries
    }

    state.counters["alphabet_size"] = seqan3::alphabet_size<alphabet_t>;
    state.counters["text_size"] = text_size;
    state.counters["k"] = k;
    state.counters["query_size"] = 150;
    state.counters["used_hashtable"] = false;
    state.counters["paralell"] = 3;
    state.counters["memory_used_in_bytes"] = sizeof(index);
}


// [SEQ] FM
template<seqan3::alphabet alphabet_t>
static void fm_search_vs_text_size(benchmark::State& state, size_t text_size)
{
    // generate text and queries
    auto gen = input_generator<alphabet_t>(1234);

    auto queries = generate_query_set<alphabet_t>(1e5);
    std::vector<alphabet_t> text = gen.generate_text(text_size, queries);

    auto index = seqan3::fm_index(text);

    size_t i = 0;
    for (auto _ : state)
    {
        benchmark::DoNotOptimize(seqan3::search(queries.at(i), index));
        i = i < queries.size()-1 ? i + 1 : 0;
    }

    state.counters["alphabet_size"] = seqan3::alphabet_size<alphabet_t>;
    state.counters["text_size"] = text_size;
    state.counters["query_size"] = 150;
    state.counters["paralell"] = 1;
    state.counters["memory_used_in_bytes"] = sizeof(index);
}


void register_all(int argc, char** argv)
{
    size_t base_size = 100000;
    for (size_t i = 1; i < 20; i++)
    {
        benchmark::RegisterBenchmark("multi_10_11_12_kmer_search", &kmer_search_vs_text_size<seqan3::dna4, 10>, i*base_size);
        benchmark::RegisterBenchmark("fm_search", &fm_search_vs_text_size<seqan3::dna4>, i*base_size);
    }

    benchmark::Initialize(&argc, argv);

    try
    {
        benchmark::RunSpecifiedBenchmarks();
    }
    catch (std::bad_alloc)
    {
        std::cerr << "out of memory. aborting safely...";
        return;
    }
}

/* EXECUTABLE ARGUMENTS:
--benchmark_format=console
--benchmark_counters_tabular=true
--benchmark_out=/home/clem/Documents/Workspace/kmer_index/source/benchmarks/fm_kmer_exact_vs_text_length/raw.csv
--benchmark_out_format=csv
--benchmark_repetitions=10
--benchmark_report_aggregates_only=true
 */

int main(int argc, char** argv)
{
    register_all(argc, argv);
    cleanup_csv("/home/clem/Documents/Workspace/kmer_index/source/benchmarks/fm_kmer_exact_vs_text_length/raw.csv");
}














// [1] https://www.illumina.com/science/technology/next-generation-sequencing/plan-experiments/read-length.html