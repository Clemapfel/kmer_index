//
// Created by Clemens Cords on 2/7/20.
// Copyright (c) 2020 Clemens Cords. All rights reserved.
//

#include <kmer_index.hpp>
#include <benchmarks/input_generator.hpp>
#include <fast_pow.hpp>

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
#include <choose_best_k.hpp>
#include <thread_pool.hpp>

int main()
{
    using namespace seqan3;

    input_generator<seqan3::dna4> input;
    auto query = input.generate_sequence(19);
    auto text = input.generate_text(1e5, {query});

    auto kmer = kmer::make_kmer_index<4, 15>(text);
    auto fm = fm_index(text);


    auto kmer_res = kmer.experimental_search(query);
    auto fm_res = seqan3::search(query, fm);



    seqan3::debug_stream << "KMER: " << kmer_res.to_vector() << "\n";
    seqan3::debug_stream << "FM  : ";
    for (auto res : fm_res)
    {
        seqan3::debug_stream << res.reference_begin_pos() << ",";
    }
    seqan3::debug_stream << "\n";
}