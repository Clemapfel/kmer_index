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
    auto text = "ACGTACGTACGTACG"_dna4;

    auto slow_index = make_slow_kmer_index<4>(text);
    auto fm = fm_index(text);

    auto query = "C"_dna4;
    debug_stream << slow_index.search(query) << "\n";
    debug_stream << search(query, fm) << "\n";

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
