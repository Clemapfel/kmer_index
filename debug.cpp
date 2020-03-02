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
    auto text = input_generator<dna4>::generate_text(1000, {"ACGT"_dna4});

    auto index = make_kmer_index<4>(text);
    fm_index fm{text};

        debug_stream << index.search("ACGT"_dna4) << " | ";
        debug_stream << search("ACGT"_dna4, fm) << "\n";
}
