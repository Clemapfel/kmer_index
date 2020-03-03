//
// Created by Clemens Cords on 3/3/20.
// Copyright (c) 2020 Clemens Cords. All rights reserved.
//

#pragma once

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


// benchmark arguments struct for readability
struct benchmark_arguments
{
    const size_t _query_size;               // size of every query
    const size_t _n_queries;                // number of queries generated
    const size_t _text_size;                // length of text

    // ctor
    benchmark_arguments(size_t query_size, size_t n_queries, size_t text_size);

    // generate queries and text pseudo-randomly
    template<seqan3::alphabet alphabet_t>
    void generate_queries_and_text(
        std::vector<std::vector<alphabet_t>>* queries,  // out
        std::vector<alphabet_t>* text,                  // out
        bool reset_state = true);                       // if true, returns same queries and text each call

    // add custom counters to keep track of benchmark arguments
    template<seqan3::alphabet alphabet_t, bool use_hashtable, int k>
    void add_counters_to(benchmark::State& state) const;
};

// #####################################################################################################################

// benchmark: kmer index exact search
template<seqan3::alphabet alphabet_t, bool use_da, size_t k>
static void kmer_search(benchmark::State& state, benchmark_arguments input);

// benchmark:: multi kmer index exact search
template<seqan3::alphabet alphabet_t, bool use_da, size_t... ks>
static void multi_kmer_search(benchmark::State& state, benchmark_arguments input);

// benchmark: fm index exact search
template<seqan3::alphabet alphabet_t>
static void fm_search(benchmark::State& state, benchmark_arguments input);

// benchmark: kmer_index construction
template<seqan3::alphabet alphabet_t, bool use_da, size_t k>
static void kmer_construction(benchmark::State& state, benchmark_arguments input);

// benchmark: multi_kmer construction
template<seqan3::alphabet alphabet_t, bool use_da, size_t... ks>
static void multi_kmer_search(benchmark::State& state, benchmark_arguments input);

// benchmark: fm_index construction
template<seqan3::alphabet alphabet_t>
static void fm_construction(benchmark::State& state, benchmark_arguments input);

//#####################################################################################################################

// function that can be expanded in fold expression
template<seqan3::alphabet alphabet_t, bool use_hashtable, size_t k>
void register_kmer_benchmarks(benchmark_arguments config);

// wrapper that programmatically registers benchmarks for all permutations of template params and args
template<seqan3::alphabet alphabet_t, bool use_hashtable, size_t... ks>
void register_all_benchmarks(std::vector<size_t> text_sizes);

// cleanup resulting csv for later parsing
void cleanup_csv(std::string path);