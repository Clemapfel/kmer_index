//
// Created by Clemens Cords on 2/7/20.
// Copyright (c) 2020 Clemens Cords. All rights reserved.
//

#include <kmer_index.hpp>

#include <seqan3/core/debug_stream.hpp>
#include <seqan3/alphabet/all.hpp>
#include <seqan3/alphabet/cigar/all.hpp>
#include <seqan3/search/fm_index/fm_index.hpp>
#include <seqan3/search/algorithm/search.hpp>

using namespace seqan3;

int main()
{
    auto text = "ACTGACTGCATCGATCGATCGTACGTACG"_dna4;
    auto single = make_kmer_index<4>(text);
    //auto collection =  make_kmer_index<dna4, 3, 4, 5>(text);
    fm_index fm{text};

    auto queries = std::vector{"ACGT"_dna4, "ACTGAC"_dna4, "GACTGCATC"_dna4, "CTGCATCGATCGATCGTACGTACG"_dna4};

    for (auto q : queries)
    {
        debug_stream << single.search(q) << " | ";
        debug_stream << search(q, fm) << "\n";
    }
}
