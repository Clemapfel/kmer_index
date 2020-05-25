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
constexpr size_t k = 5;
constexpr size_t text_size = 10000;

int main()
{
    auto query = "ACGTAACGTA"_dna4;

    auto input = input_generator<dna4>(0);
    auto text = input.generate_text(text_size, {query});
    auto kmer = kmer_index_element<seqan3::dna4, k, uint32_t>(text);

    auto result_true = kmer.search(query).to_vector();
    auto result_test = kmer.search_test(query).to_vector();

    seqan3::debug_stream << result_true << "\n\n" << result_test << "\n\n";

    auto fm = fm_index(text);
    std::vector<size_t> fm_result;
    for (auto _ : search(query, fm))
        fm_result.push_back(_.second);

    seqan3::debug_stream << fm_result << "\n";

    seqan3::debug_stream << "kmer_old : " << result_true.size() << "\nkmer_new : " << result_test.size() << "\nfm : " << fm_result.size();
}

