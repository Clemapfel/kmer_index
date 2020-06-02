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

#include <benchmarks/benchmarks.hpp>

using namespace seqan3;
using alphabet_t = dna15;
constexpr size_t k = 10;
constexpr size_t text_size = 100000;

void force_error(std::vector<alphabet_t> query, size_t seed)
{
    auto input = input_generator<alphabet_t>(seed);
    auto text = input.generate_text(text_size, {});
    auto kmer = kmer::make_kmer_index<k>(text);
    auto fm = fm_index(text);

    std::vector<unsigned int> fm_result;
    for (auto _ : search(query, fm))
        fm_result.push_back(_.second);

    auto res = kmer.search<k>(query);
    auto kmer_result = res.to_vector();

    auto equal = fm_result == kmer_result;

    seqan3::debug_stream << (equal ? "EQUAL FOR " : "NOT EQUAL FOR ") << "\nQUERY " << query << " (" << query.size() << ")\n"
                         << "\nFM : " << fm_result << "\n\n KMER : " << kmer_result << "\n";
    seqan3::debug_stream << "query size = " << query.size() << "\nseed = " << seed << "\n"
                         <<  "difference (fm - kmer) = " << int(fm_result.size()) - int(kmer_result.size()) << "\n";

    std::vector<uint32_t> diff;
    std::set_difference(fm_result.begin(), fm_result.end(), kmer_result.begin(), kmer_result.end(), std::inserter(diff, diff.begin()));
    seqan3::debug_stream << diff;

    if (not equal)
        exit(1);
    else
        exit(0);
}

int main()
{
    benchmark_input<seqan3::dna4> input(1000, 1000);

    benchmark::RegisterBenchmark("test", &kmer_search<seqan3::dna4, 10>);//, input);

    force_error("Y"_dna15, 5001);
    return 0;

    /*
    for (size_t i = 0; i < 10000; ++i)
    {
        auto input = input_generator<dna4>(i);
        auto text = input.generate_text(text_size, {});
        auto kmer = kmer_index_element<seqan3::dna4, k, uint32_t>(text);
        auto fm = fm_index(text);

        for (size_t query_size : {1ul, k-2, k-1, k, k+1, k+2, 2*k, 2*k+1, 3*k, 5*k, 10*k})
        {
            auto query = input.generate_sequence(query_size);

            std::vector<unsigned int> fm_result;
            for (auto _ : search(query, fm))
                fm_result.push_back(_.second);

            auto kmer_result = kmer.search(query).to_vector();

            auto equal = fm_result == kmer_result;

            if (equal)
                seqan3::debug_stream << "TRUE ";
            else
            {
                seqan3::debug_stream << "\n\nNOT EQUAL FOR QUERY " << query << " (" << query.size() << ")\n"
                                     << "\nFM : " << fm_result << "\n\n KMER : " << kmer_result << "\n";
                seqan3::debug_stream << "query size = " << query_size << "\nseed = " << i << "\n"
                <<  "difference (fm - kmer) = " << int(fm_result.size()) - int(kmer_result.size()) << "\n";

                std::vector<uint32_t> diff;
                std::set_difference(fm_result.begin(), fm_result.end(), kmer_result.begin(), kmer_result.end(), std::inserter(diff, diff.begin()));
                seqan3::debug_stream << diff;
                exit(1);
            }
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

