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


    int main()
    {
        auto query = "ACGT"_dna4;
        auto input = input_generator<dna4>(1235);
        std::vector<std::vector<dna4>> queries = input.generate_queries(100, 6);

        auto text = input.generate_sequence(1e6);
        auto map_index = make_kmer_index<5, 6, 7>(text);

        debug_stream << map_index.calculate_size();
        //debug_stream << map_index.search(query).size();
        return 0;


        /*
        debug_stream << "starting test...\n";

        for (size_t i = 0; i < 1; ++i)
        {
            // state not reset so new text everytime
            auto text = input.generate_sequence(1e3);
            auto da_index = make_kmer_index<true, 5, 6, 7>(text);
            auto map_index = make_kmer_index<false, 5, 6, 7>(text);
            auto fm = fm_index(text);

            bool results_equal;
            try
            {
                auto fm_results = search(query, fm);
                size_t fm_size = 0;
                for (auto f : fm_results)
                    fm_size++;

                auto size = da_index.search(query).size() + map_index.search(query).size() + fm_size;
                //auto fm_results = search(query, fm);
                size = da_index.search(query).size() + map_index.search(query).size() + fm_size;
                results_equal = (size - (fm_size * 3)) == 0;
            }
            catch (std::out_of_range)
            {
                results_equal = false;
                std::cerr << "out of range exception\n";
            }

            if (not results_equal)
            {
                debug_stream << "results not equal for i = " << i << "\n";
                //debug_stream << text << "\n\n";
                debug_stream << "fm  : " << search(query, fm) << "\n";
                debug_stream << "da  : " << da_index.search(query) << "\n";
                debug_stream << "map : " << map_index.search(query) << "\n";

                return 1;
            }

        }

        debug_stream << "test passed succesfully.";
         */
    }




