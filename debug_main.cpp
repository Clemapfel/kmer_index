//
// Created by Clemens Cords on 2/7/20.
// Copyright (c) 2020 Clemens Cords. All rights reserved.
//

#include <kmer_index.hpp>
#include <input_generator.hpp>

#include <seqan3/core/debug_stream.hpp>
#include <seqan3/alphabet/all.hpp>
#include <seqan3/alphabet/cigar/all.hpp>
#include <seqan3/search/fm_index/fm_index.hpp>
#include <seqan3/search/algorithm/search.hpp>

#include <thread>
#include <chrono>
#include <memory>

using namespace seqan3;

int main()
{
    std::vector<std::vector<dna4>> queries = input_generator<dna4, 1234>::generate_queries(100, 6);
    debug_stream << "starting test...\n";

    for (size_t i = 0; i < 1; ++i)
    {
        // state not reset so new text everytime
        auto text = input_generator<dna4, 1234>::generate_sequence(1e3);

        auto da_index = make_kmer_index<true, 5, 6, 7>(text);
        //auto map_index = make_kmer_index<false, 5, 6, 7>(text);
        //auto fm = fm_index(text);

        // search in paralell
        auto paralell_results = da_index.search_par(queries);

        debug_stream << "paralell done.\n";

        // search in sequence
        std::vector<std::vector<unsigned int>> sequence_results;
        size_t in = 0;
        for (auto& q : queries) {
            debug_stream << "starting serach " << std::to_string(in++) << "\n";
            da_index.search_seq(q);

            TODO: why does it work in paralell but not in sequence?

        }

        debug_stream << "seq done.\n";

        std::cout << (paralell_results == sequence_results ? "test passed" : "test failed");

        /*
        bool results_equal;
        try
        {
            //auto fm_results = search(query, fm);
             auto size = da_index.search(query).size() + map_index.search(query).size() + fm_results.size();
             results_equal = (size - (fm_results.size() * 3)) == 0;
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
         */
    }

    debug_stream << "test passed succesfully.";
}
