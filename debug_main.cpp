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
    size_t _min_hash = 0, _max_hash = 1023;
    _shift_amount = 64 - log2l(_max_hash - _min_hash);
    seqan3::debug_stream << "shift amount: " << _shift_amount << "\n";

    size_t min = std::numeric_limits<size_t>::max(), max = 0;
    std::vector<size_t> hash_hashes;

    for (size_t i = _min_hash; i < _max_hash; ++i)
    {
        auto h = hash_hash(i);
        hash_hashes.push_back(h);

        if (h < min)
            min = h;

        if (h > max)
            max = h;
    }

    std::sort(hash_hashes.begin(), hash_hashes.end());
    seqan3::debug_stream << hash_hashes;
    seqan3::debug_stream << "\n_____________________ \n";
    seqan3::debug_stream << "n: " << hash_hashes.begin() << " | min: " << min << " | max: " << max << "\n";


    /*
    auto query = "ACGT"_dna4;

    std::cout << "starting test...\n";

    for (size_t i = 0; i < 1; ++i)
    {
        // state not reset so new text everytime
        auto text = input_generator<dna4, 1234>::generate_sequence(1000);

        auto index = make_kmer_index<true, 5>(text);
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
