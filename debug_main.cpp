//
// Created by Clemens Cords on 2/7/20.
// Copyright (c) 2020 Clemens Cords. All rights reserved.
//

#include <kmer_index.hpp>
#include <benchmarks/input_generator.hpp>

#include <seqan3/core/debug_stream.hpp>
#include <seqan3/alphabet/all.hpp>
#include <seqan3/alphabet/cigar/all.hpp>
#include <seqan3/search/fm_index/fm_index.hpp>
#include <seqan3/search/search.hpp>

#include <thread>
#include <chrono>
#include <memory>
#include <future>
#include <functional>
#include <thread_pool.hpp>

using namespace seqan3;
using alphabet_t = dna4;
constexpr size_t k = 6;
constexpr size_t text_size = 100000;

int main()
{
    auto input = input_generator<dna4>(0);
    auto text = input.generate_sequence(text_size);
    auto kmer = kmer_index_element<seqan3::dna4, k, uint32_t>(text);

    auto result = kmer.search("AGCT"_dna4);

    seqan3::debug_stream << result.to_vector() << "\n";

    debug_stream << "\nstarting loop:\n";
    size_t previous = -1;   //sic
    for (auto pos : result)
    {
        if (previous == pos)
            return 1;

        debug_stream << pos << " | ";
        previous = pos;
    }

    debug_stream << "\n";
}

