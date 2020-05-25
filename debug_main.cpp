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
    for (size_t i = 0; i < 100; ++i)
    {
        auto input = input_generator<dna4>(i);
        auto text = input.generate_text(text_size, {});
        auto kmer = kmer_index_element<seqan3::dna4, k, uint32_t>(text);

        auto query = input.generate_sequence(2*k);

        auto fm = fm_index(text);
        std::vector<unsigned int> fm_result;
        for (auto _ : search(query, fm))
            fm_result.push_back(_.second);

        auto kmer_result = kmer.search_test(query).to_vector();

        auto equal = fm_result == kmer_result;

        if (equal)
            seqan3::debug_stream << "TRUE ";
        else
        {
            seqan3::debug_stream << "NOT EQUAL FOR QUERY " << query
            << "\nFM : " << fm_result << "\n\n KMER : " << kmer_result << "\n";
            exit(1);
        }
    }

    return 0;

    /*
    auto query = "ACGTAACGTA"_dna4;



    auto result_test = kmer.search_test(query).to_vector();
    //auto result_true = kmer.search(query).to_vector();

    //seqan3::debug_stream << result_true << "\n\n";
    seqan3::debug_stream << result_test << "\n\n";

    //query = "ACGTAACGTAAC"_dna4;

    auto fm = fm_index(text);
    std::vector<size_t> fm_result;
    for (auto _ : search(query, fm))
        fm_result.push_back(_.second);

    seqan3::debug_stream << fm_result << "\n";
    //seqan3::debug_stream << result_test << "\n";

    //seqan3::debug_stream << "kmer_old : " << result_true.size() << "\nkmer_new : " << result_test.size() << "\nfm : " << fm_result.size();
    */
}

