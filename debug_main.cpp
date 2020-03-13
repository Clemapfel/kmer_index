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

uint8_t _shift_amount;
size_t hash_hash(size_t hash)
{
    return (hash * 11400714819323198485llu) >> _shift_amount;
}

int main()
{
    /*
    size_t _min_hash = 0, _max_hash = 1023;
    _shift_amount = 64 - log2l(_max_hash - _min_hash) + 1;
    seqan3::debug_stream << "shift amount: " << _shift_amount << "\n";

    size_t min = std::numeric_limits<size_t>::max(), max = 0;
    std::set<size_t> hash_hashes;

    for (size_t i = _min_hash; i < _max_hash; ++i)
    {
        auto h = hash_hash(i);
        hash_hashes.insert(h);

        if (h < min)
            min = h;

        if (h > max)
            max = h;
    }

    seqan3::debug_stream << hash_hashes;
    seqan3::debug_stream << "\n_____________________ \n";
    seqan3::debug_stream << "n: " << hash_hashes.size() << " | min: " << min << " | max: " << max << "\n";
     */
    auto query = "ACGT"_dna4;

    std::cout << "starting test...\n";

    for (size_t i = 0; i < 1; ++i)
    {
        // state not reset so new text everytime
        auto text = input_generator<dna4, 1234>::generate_sequence(1e6);

        auto da_index = make_kmer_index<true, 8>(text);
        auto map_index = make_kmer_index<false, 8>(text);
        auto fm = fm_index(text);

        bool results_equal;
        try
        {
            auto fm_results = search(query, fm);
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
    }

    std::cout << "test passed succesfully.";
}
