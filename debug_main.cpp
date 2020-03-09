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

using namespace seqan3;

int main()
{
    auto query = "ACGT"_dna4;

    std::cout << "starting test...\n";

    for (size_t i = 0; i < 1000; ++i)
    {
        // state not reset so new text everytime
        auto text = input_generator<dna4, 1234>::generate_sequence(1e3);

        auto index = make_fast_kmer_index<3>(text);
        auto fm = fm_index(text);

        bool results_equal;
        try
        {
             results_equal = (index.search(query).size() == search(query, fm).size());
        }
        catch (std::out_of_range)
        {
            results_equal = false;
            std::cerr << "out of range exception\n";
        }

        if (not results_equal)
        {
            debug_stream << "results not equal for i = " << i << "\n";
            debug_stream << text << "\n\n";
            debug_stream << "fm  : " << search(query, fm) << "\n";
            debug_stream << "kmer: " << index.search(query) << "\n";
            return 1;
        }
    }

    std::cout << "test passed succesfully.";

    /*
    for (auto i : std::vector{15, 14, 13, 12, 11, 10, 9, 8, 7, 6, 5, 4})
    {
        auto query = input_generator<dna4>::generate_sequence(i);

        debug_stream << "<" << i << ">\n";
        debug_stream << slow_index.search(query) << "\n";
        debug_stream << search(query, fm) << "\n\n";
    }
     */
}
