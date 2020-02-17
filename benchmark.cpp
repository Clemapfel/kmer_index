//
// Created by Clemens Cords on 2/7/20.
// Copyright (c) 2020 Clemens Cords. All rights reserved.
//
#include <vector>
#include <type_traits>
#include <sstream>

#include <benchmark/benchmark.h>

#include <seqan3/alphabet/concept.hpp>
#include <seqan3/alphabet/aminoacid/all.hpp>
#include <seqan3/alphabet/nucleotide/all.hpp>
#include <seqan3/search/algorithm/search.hpp>
#include <seqan3/search/fm_index/fm_index.hpp>
#include <seqan3/std/type_traits>

#include "input_generator.hpp"
#include "kmer_index.hpp"

// convert to string for naming
template<seqan3::alphabet T>
std::string alphabet_to_string()
{
    if (std::is_same<T, seqan3::dna4>::value)
        return "dna4";
    else if (std::is_same<T, seqan3::dna5>::value)
        return "dna5";
    else if (std::is_same<T, seqan3::dna15>::value)
        return "dna15";
    else if (std::is_same<T, seqan3::aa10li>::value or std::is_same<T, seqan3::aa10murphy>::value)
        return "aa10";
    else if (std::is_same<T, seqan3::aa20>::value)
        return "aa20";
    else if (std::is_same<T, seqan3::aa27>::value)
        return "aa27";
    else if (std::is_same<T, seqan3::rna4>::value)
        return "rna4";
    else if (std::is_same<T, seqan3::rna5>::value)
        return "rna5";
    else if (std::is_same<T, seqan3::rna15>::value)
        return "rna15";
    else
        return "other";
}

// benchmark arguments struct for readability
struct benchmark_arguments
{
    const size_t _query_size;               // size of every query
    const size_t _n_queries;                // number of queries generated
    const size_t _text_size;                // length of text
    const float _no_hit_queries_ratio;      // text containts ratio * n_queries of the queries, 0 hits for rest

    // ctor
    benchmark_arguments(
        size_t query_size,
        size_t n_queries,
        size_t text_size,
        float no_hit_queries_ratio = 0)
        : _query_size(query_size),
          _n_queries(n_queries),
          _text_size(text_size),
          _no_hit_queries_ratio(no_hit_queries_ratio)
    {
    }

    // generate queries and text pseudo-randomly
    template<seqan3::alphabet alphabet_t>
    void generate_queries_and_text(
        std::vector<std::vector<alphabet_t>>* queries,  // out
        std::vector<alphabet_t>* text,                  // out
        bool reset_state = true)                        // if true, returns same queries and text each call
    {
        if (reset_state)
            input_generator<alphabet_t>::reset_state();

        queries->clear();
        for (auto q : input_generator<alphabet_t>::generate_queries(_n_queries, _query_size))
            queries->push_back(q);

        text->clear();
        for (auto c : input_generator<alphabet_t>::generate_text(_text_size, *queries))
            text->push_back(c);
    }

    template<seqan3::alphabet alphabet_t, bool use_da, int k>
    void add_counters_to(benchmark::State& state) const
    {
        state.counters["alphabet_size"] = seqan3::alphabet_size<alphabet_t>;
        state.counters["text_size"] = _text_size;
        state.counters["k"] = k;
        state.counters["query_size"] = _query_size;
        state.counters["n_queries"] = _n_queries;
        state.counters["no_hit_query_ratio"] = _no_hit_queries_ratio;
        state.counters["used_da"] = use_da;
    }

    // generate csv header
    static const std::string get_header()
    {
        return "name,alphabet,text_size,k,query_size,n_queries,no_hit_query_ratio,uses_direct_addressing,report_type";
    }

    // generate benchmark name which forms first few columns of csv
    const std::string get_name(std::string id, std::string alphabet_id, bool use_da, int k) const
    {
        return
           "" + id + ","
           + "" + alphabet_id + ","
           + std::to_string(_text_size) + ","
           + (k != 0 ? std::to_string(k) : "NaN") + ","
           + std::to_string(_query_size) + ","
           + std::to_string(_n_queries )+ ","
           + std::to_string(_no_hit_queries_ratio) + ","
           + std::to_string(use_da) + ",";
    }
};

// #####################################################################################################################

// benchmark: kmer index exact search
template<seqan3::alphabet alphabet_t, bool use_da, size_t k>
static void kmer_search(benchmark::State& state, benchmark_arguments input)
{
    // generate text and queries
    std::vector<std::vector<alphabet_t>> queries;
    std::vector<alphabet_t> text;
    input.generate_queries_and_text(&queries, &text, true);

    kmer_index<alphabet_t, k, uint64_t, uint32_t, use_da> index{text};

    // add custom vars to output
    input.add_counters_to<alphabet_t, use_da, k>(state);
    state.counters["memory_used(mb)"] = sizeof(index) / 1e6;

    // cycle through queries
    unsigned long long i = 0;
    for (auto _ : state)
    {
        benchmark::DoNotOptimize(index.search(queries.at(i)));
        i = i < queries.size()-1 ? i + 1 : 0;
    }
}

/*
// benchmark:: multi kmer index exact search
template<seqan3::alphabet alphabet_t, bool use_da, size_t... ks>
static void multi_kmer_search(benchmark::State& state, benchmark_arguments input)
{
    std::vector<std::vector<alphabet_t>> queries;
    std::vector<alphabet_t> text;
    input.generate_queries_and_text(&queries, &text, true);

    multi_kmer_index<alphabet_t, k, uint64_t, uint32_t, use_da, ks...> index{text};

    unsigned long long i = 0;
    for (auto _ : state)
    {
        benchmark::DoNotOptimize(index.search(queries.at(i)));
        i = i < queries.size()-1 ? i + 1 : 0;
    }
}
*/

// benchmark: fm index exact search
template<seqan3::alphabet alphabet_t>
static void fm_search(benchmark::State& state, benchmark_arguments input)
{
    std::vector<std::vector<alphabet_t>> queries;
    std::vector<alphabet_t> text;
    input.generate_queries_and_text(&queries, &text, true);

    seqan3::fm_index index{text};

    input.add_counters_to<alphabet_t, false, -1>(state);
    state.counters["memory_used(mb)"] = sizeof(index) / 1e6;

    unsigned long long i = 0;
    for (auto _ : state)
    {
        benchmark::DoNotOptimize(search(queries.at(i), index));
        i = i < queries.size()-1 ? i + 1 : 0;
    }
}

// benchmark: kmer_index construction
template<seqan3::alphabet alphabet_t, bool use_da, size_t k>
static void kmer_construction(benchmark::State& state, benchmark_arguments input)
{
    std::vector<std::vector<alphabet_t>> queries;
    std::vector<alphabet_t> text;
    input.generate_queries_and_text<seqan3::dna4>(&queries, &text, true);

    for (auto _ : state)
        benchmark::DoNotOptimize(kmer_index<alphabet_t, k, uint64_t, uint32_t, use_da>{text});

    // log memory used
    kmer_index<alphabet_t, k, uint64_t, uint32_t, use_da> index{text};

    input.add_counters_to<alphabet_t, use_da, k>(state);
    state.counters["memory_used(mb)"] = sizeof(index) / 1e6;
}

/*
// benchmark: multi_kmer construction
template<seqan3::alphabet alphabet_t, bool use_da, size_t... ks>
static void multi_kmer_search(benchmark::State& state, benchmark_arguments input)
{
    std::vector<std::vector<alphabet_t>> queries;
    std::vector<alphabet_t> text;-
    input.generate_queries_and_text(&queries, &text, true);

    for (auto _ : state)
        benchmark::DoNotOptimize(multi_kmer_index<alphabet_t, k, uint64_t, uint32_t, use_da, ks...> index{text});
}
 */

// benchmark: fm_index construction
template<seqan3::alphabet alphabet_t>
static void fm_construction(benchmark::State& state, benchmark_arguments input)
{
    std::vector<std::vector<alphabet_t>> queries;
    std::vector<alphabet_t> text;
    input.generate_queries_and_text<seqan3::dna4>(&queries, &text, true);

    for (auto _ : state)
        benchmark::DoNotOptimize(seqan3::fm_index{text});

    // log memory used
    seqan3::fm_index index{text};
    input.add_counters_to<alphabet_t, false, -1>(state);
    state.counters["memory_used(mb)"] = sizeof(index) / 1e6;
}

//#####################################################################################################################

// wrapper that programmatically registers benchmarks
template<seqan3::alphabet alphabet_t, bool use_da, size_t... ks>
class register_benchmarks
{
    public:
        // initialize configs with all permutations of:
        static void register_all(
            std::map<size_t, std::vector<size_t>> query_sizes_per_k,
            std::vector<size_t> text_sizes,
            std::vector<float> no_hit_ratios)
        {
            _configs.clear();

            for (auto k : std::vector{ks...})
            {
                std::vector<size_t> query_sizes;
                if (query_sizes_per_k.find(k) != query_sizes_per_k.end())
                    query_sizes = query_sizes_per_k.at(k);
                else
                    query_sizes = {k};

                for (auto query_size : query_sizes)
                {
                    for (auto text_size : text_sizes)
                    {
                        for (auto no_hit_ratio : no_hit_ratios)
                        {
                            size_t n_queries = 0.1 * text_size > 100 ? 0.1 * text_size : 100;

                            if (_configs.find(k) == _configs.end())
                                _configs.insert(std::make_pair(k, std::vector<benchmark_arguments>{}));

                            _configs[k].emplace_back(query_size, n_queries, text_size, no_hit_ratio);
                        }
                    }
                }
            }

            for (auto pair : _configs)
            {
                for (auto& config : pair.second)
                {
                    // add replaced later on
                    std::string name = "kmer_exact_search";
                    (benchmark::RegisterBenchmark(name.c_str(), &kmer_search<alphabet_t, use_da, ks>, config), ...);

                    name = "fm_exact_search";
                    benchmark::RegisterBenchmark(name.c_str(), &fm_search<alphabet_t>, config);

                    name = "kmer_construction";
                    (benchmark::RegisterBenchmark(name.c_str(), &kmer_construction<alphabet_t, use_da, ks>, config), ...);

                    name = "fm_construction";
                    benchmark::RegisterBenchmark(name.c_str(), &fm_construction<alphabet_t>, config);
                }
            }
        }

    private:
        // fields
        static inline std::map<size_t, std::vector<benchmark_arguments>> _configs{};
        static inline std::string _alphabet_id = alphabet_to_string<alphabet_t>();
};

void cleanup_csv(std::string path)
{
    std::ifstream file_in;
    file_in.open(path);

    if (file_in.fail())
        std::cerr << "Error opening file " + path + "\naborting...";

    // name outfile
    auto it = path.find("_raw");
    std::string out_path = path.erase(it, 4);

    std::ofstream file_out;
    file_out.open(out_path);

    if (file_out.fail())
        std::cerr << "Error opening file " + out_path + "\naborting...";

    std::string line;
    bool past_header = false;

    while (std::getline(file_in, line))
    {
        if (not past_header)
        {
            // csv header line
            if (line.substr(0, 4) == "name")
            {
                line.erase(std::remove(line.begin(), line.end(), '\"'), line.end());
                file_out << line;
                past_header = true;
            }
            // skip gbenchmark header
            else
                continue;
        }
        else
            file_out << line;
    }

    std::cout << "[DEBUG] csv written to " + out_path + "\n";

    file_in.close();
    file_out.close();
}

// #####################################################################################################################

// main
int main(int argc, char** argv)
{
    std::map<size_t, std::vector<size_t>> query_sizes_per_k{};
    query_sizes_per_k[4] = {3, 4, 5};
    query_sizes_per_k[5] = {5};
    query_sizes_per_k[6] = {5, 6};

    register_benchmarks<seqan3::dna4, false, 4, 5, 6>::register_all(
            query_sizes_per_k,
            // text size
            {1000, 3000},
            // hit ratios
            {1}
    );

    register_benchmarks<seqan3::dna4, true, 5>::register_all(
            {{5, {5}}},
            {3000},
            {1}
    );

    benchmark::Initialize(&argc, argv);
    benchmark::RunSpecifiedBenchmarks();

    cleanup_csv("../source/benchmark_out_raw.csv");
}


