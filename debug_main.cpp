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
#include <benchmarks/input_generator.hpp>
#include <thread_pool.hpp>
#include <kmer_index_result.hpp>

using namespace seqan3;


constexpr size_t k = 7;
int main()
{
    const auto results = std::vector<uint32_t>{4, 5, 123, 123, 154, 167, 132};

    auto not_flat = std::vector<const std::vector<uint32_t>*>{&results, &results, &results};

    kmer_index_results<dna4, k> results_flat{&results};
    kmer_index_results<dna4, k> results_not_flat{not_flat};

    debug_stream << results_flat.to_vector() << "\n";
    debug_stream << results_not_flat.to_vector() << "\n";




    /*
    input_generator<seqan3::dna4> input;
    auto text = input.generate_sequence(100000);

    auto kmer = kmer_index<seqan3::dna4, k>{text};
    auto fm = fm_index{text};

    std::vector<std::vector<seqan3::dna4>> k_queries = input.generate_queries(1000, k);
    auto subk_queries = input.generate_queries(1000, k-2);
    auto nk_queries = input.generate_queries(1000, 4*k);

    for (auto q : k_queries)
    {
        auto kmer_results = kmer.search(q);
        auto fm_results = seqan3::search(q, fm);

        std::sort(kmer_results.begin(), kmer_results.end());

        size_t fm_size = 0;
        for (auto _ : fm_results)
            fm_size += 1;

        if (not kmer_results.size() == fm_size)
            return 1;
    }

    return 0;
     */
}




