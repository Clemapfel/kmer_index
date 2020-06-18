#include <kmer_index.hpp>
#include <benchmarks/input_generator.hpp>
#include <fast_pow.hpp>
#include <benchmarks/cleanup_csv.hpp>

#include <benchmark/benchmark.h>

#include <seqan3/alphabet/nucleotide/dna4.hpp>
#include <seqan3/search/fm_index/fm_index.hpp>
#include <seqan3/search/fm_index/bi_fm_index.hpp>
#include <seqan3/search/search.hpp>

#include <numeric>
#include <random>

size_t start_seed = 1234;
size_t element_seed = start_seed;
size_t single_seed = start_seed;
size_t multi_seed = start_seed;

// [SEQ] single kmer search
template<seqan3::alphabet alphabet_t, size_t k>
static void element(benchmark::State& state, size_t text_length, size_t query_length)
{
    state.counters["seed"] = element_seed;

    input_generator<alphabet_t> input(element_seed);
    auto text = input.generate_sequence(text_length);
    auto queries = input.generate_queries(200000, query_length);

    auto index = kmer::detail::kmer_index_element<alphabet_t, unsigned int, k>(text);

    size_t i = 0;
    for (auto _ : state)
    {
        auto res = index.search(queries.at(i));
        i = i < queries.size() - 1 ? i + 1 : 0;
    }

    state.counters["text_length"] = text_length;
    state.counters["query_length"] = query_length;
    state.counters["alphabet_size"] = seqan3::alphabet_size<alphabet_t>;
    state.counters["k_0"] = k;
    element_seed++;
}

// [SEQ] multi kmer search
template<seqan3::alphabet alphabet_t, size_t k>
static void multi_single(benchmark::State& state, size_t text_length, size_t query_length)
{
    state.counters["seed"] = single_seed;

    input_generator<alphabet_t> input(single_seed);
    auto text = input.generate_sequence(text_length);
    auto queries = input.generate_queries(200000, query_length);

    auto index = kmer::make_kmer_index<k>(text);

    size_t i = 0;
    for (auto _ : state)
    {
        auto res = index.kmer::detail::template kmer_index_element<alphabet_t, uint32_t, k>::search(queries.at(i));
        i = i < queries.size() - 1 ? i + 1 : 0;
    }

    state.counters["text_length"] = text_length;
    state.counters["query_length"] = query_length;
    state.counters["alphabet_size"] = seqan3::alphabet_size<alphabet_t>;
    state.counters["k_0"] = k;
    single_seed++;
}

// [SEQ] multi kmer search
template<seqan3::alphabet alphabet_t, size_t k>
static void multi_multi(benchmark::State& state, size_t text_length, size_t query_length)
{
    state.counters["seed"] = multi_seed;

    input_generator<alphabet_t> input(multi_seed);
    auto text = input.generate_sequence(text_length);
    auto queries = input.generate_queries(200000, query_length);

    auto index = kmer::make_kmer_index<k>(text);

    size_t i = 0;
    for (auto _ : state)
    {
        auto res = index.search(queries.at(i));
        i = i < queries.size() - 1 ? i + 1 : 0;
    }

    state.counters["text_length"] = text_length;
    state.counters["query_length"] = query_length;
    state.counters["alphabet_size"] = seqan3::alphabet_size<alphabet_t>;
    state.counters["k_0"] = k;
    multi_seed++;
}

int main(int argc, char** argv)
{
    size_t text_length = 300000;

    for (size_t query_length = 7; query_length <= 100; ++query_length)
    {
        benchmark::RegisterBenchmark("element", &multi_single<seqan3::dna4, 10>, text_length, query_length);
        benchmark::RegisterBenchmark("multi_single", &multi_single<seqan3::dna4, 10>, text_length, query_length);
        benchmark::RegisterBenchmark("multi_multi", &multi_multi<seqan3::dna4, 10>, text_length, query_length);
    }

    benchmark::Initialize(&argc, argv);
    benchmark::RunSpecifiedBenchmarks();

    cleanup_csv("/home/clem/Documents/Workspace/kmer_index/source/benchmarks/multi_k_overhead/raw.csv");

}