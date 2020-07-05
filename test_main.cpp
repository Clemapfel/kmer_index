//
// Created by clem on 5/29/20.
// Copyright (c) 2020 Clemens Cords. All rights reserved.
//
#include <seqan3/core/debug_stream.hpp>
#include <seqan3/alphabet/all.hpp>

#include <kmer_index.hpp>
#include <benchmarks/input_generator.hpp>
#include <seqan3/search/fm_index/fm_index.hpp>
#include <seqan3/search/search.hpp>

using alphabet_1 = seqan3::dna4;
using alphabet_2 = seqan3::dna15;

constexpr size_t k_0 = 1, k_1 = 5, k_2 = 10;
constexpr size_t text_size = 1000000;

static size_t seed = 0;

template<seqan3::alphabet alphabet_t, size_t k>
void run_test()
{
    for (size_t i = 0; i < 1000; ++i)
    {
        auto input = input_generator<alphabet_t>(seed++);
        auto text = input.generate_sequence(text_size);
        auto single_kmer = kmer::make_kmer_index<k>(text);
        auto multi_kmer = kmer::make_kmer_index<k, k+1, k+2>(text);
        auto fm = seqan3::fm_index(text);

        for (size_t query_size = k-5; query_size < 2*k; query_size++)
        {
            auto query = input.generate_sequence(query_size);

            // get results as vectors
            std::vector<unsigned int> fm_result;
            for (auto res : seqan3::search(query, fm))
                fm_result.push_back(res.reference_begin_pos());

            std::vector<unsigned int> single_kmer_result = single_kmer.search(query).to_vector();
            std::vector<unsigned int> multi_kmer_result = multi_kmer.search(query).to_vector();

            // compare
            bool equal = (fm_result == single_kmer_result) and (fm_result == multi_kmer_result);

            if (not equal)
            {
                seqan3::debug_stream << "NOT EQUAL FOR " << "\nQUERY " << query << " (" << query.size() << ")\n";
                seqan3::debug_stream << "query size = " << query.size() << "\nseed = " << seed << "\n"
                                     << "difference (fm - single) = " << int(fm_result.size()) - int(single_kmer_result.size()) << "\n"
                                     << "difference (fm - multi) = " << int(fm_result.size()) - int(multi_kmer_result.size()) << "\n";

                /*
                std::vector<uint32_t> single_diff;
                std::set_difference(fm_result.begin(), fm_result.end(), single_kmer_result.begin(), single_kmer_result.end(),
                                    std::inserter(diff, diff.begin()));
                seqan3::debug_stream << diff;*/

                exit(1);
            }
            else
                seqan3::debug_stream << "TRUE ";
        }
    }


    return;
}

// TODO: rewrite in google test
int main()
{
    seqan3::debug_stream << "starting test...\n";

    run_test<alphabet_2, k_2>();
    run_test<alphabet_2, k_1>();
    run_test<alphabet_2, k_0>();

    seqan3::debug_stream << "\ntest succesfull.\n";
}
