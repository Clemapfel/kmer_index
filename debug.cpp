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
    kmer_index<dna4, 4> single{text};
    kmer_index<dna4, 4, uint32_t, uint16_t, true> single_da{text};

    multi_kmer_index<dna4, 3, 4, 5> collection{text};
    fm_index fm{text};

    auto queries = std::vector{"ACGT"_dna4, "ACTGAC"_dna4, "GACTGCATC"_dna4, "CTGCATCGATCGATCGTACGTACG"_dna4};

    for (auto q : queries)
    {
        debug_stream << single.search(q) << " | ";
        debug_stream << collection.search(q) << " | ";
        debug_stream << search(q, fm) << "\n";
    }
}
