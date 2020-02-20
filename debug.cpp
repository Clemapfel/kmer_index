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
    std::vector<std::vector<dna4>> queries;
    for (size_t length : {21, 25, 26})
        for (auto q : input_generator<dna4>::generate_queries(1, length))
            queries.push_back(q);

    auto text = input_generator<dna4>::generate_text(10000, queries);

    auto index = make_kmer_index<5, 7>(text);
    fm_index fm{text};

    for (auto q : queries)
    {
        debug_stream << index.search(q) << " | ";
        debug_stream << search(q, fm) << "\n";
    }
}
