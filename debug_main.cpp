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
    constexpr size_t k = 5;

    void force_error(std::vector<dna4> query, size_t seed)
    {
        auto input = input_generator<dna4>(seed);
        auto text = input.generate_sequence(1e3);
        text.insert(text.begin(), query.begin(), query.end());
        auto kmer = kmer_index_element<seqan3::dna4, k, uint32_t>(text);

        auto fm = fm_index(text);
        debug_stream << "fm  : " << search(query, fm) << "\n";

        auto kmer_results = kmer.search(query);
        std::sort(kmer_results.begin(), kmer_results.end());

        debug_stream << "kmer : " << kmer_results << "\n";

        exit(0);
    }

    int main()
    {
        force_error("TAGGCCTGGAGCGTG"_dna4, 1234);

        debug_stream << "starting test...\n";

        auto input = input_generator<dna4>();

        for (size_t i = 0; i < 1000; ++i)
        {
            // state not reset so new text everytime
            input.reset_state(i);
            auto text = input.generate_sequence(1e3);
            auto kmer = kmer_index_element<seqan3::dna4, k, uint32_t>(text);
            auto fm = fm_index(text);

            for (size_t size : {k-2, k-1, k, 2*k, 3*k, 2*k+3})
            {
                //auto query = "GGCAGCATCT"_dna4;
                auto query = input.generate_sequence(size);

                auto fm_results = search(query, fm);
                size_t fm_size = 0;
                for (auto f : fm_results)
                    fm_size++;

                bool results_equal = (kmer.search(query).size() - fm_size) == 0;

                auto kmer_results =  kmer.search(query);
                std::sort(kmer_results.begin(), kmer_results.end());

                if (not results_equal)
                {

                    debug_stream << "results not equal for seed = " << i << "\n";
                    debug_stream << "text (head) : " << std::vector<dna4>(text.begin(), text.begin() + 100) << "\n";
                    debug_stream << "query : " << query << " (" << query.size() << ")\n";
                    debug_stream << "fm  : " << search(query, fm) << "\n";
                    debug_stream << "kmer : " << kmer_results << "\n";

                    return 1;
                }
            }

        }

        debug_stream << "test passed succesfully.";

    }




