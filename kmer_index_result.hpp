// Copyright (c) 2020 Clemens Cords. All rights reserved.

#pragma once

#include <kmer_index.hpp>

template<seqan3::alphabet alphabet_t, size_t k>
kmer_index_results
{


    private:
        const std::vector<position_t>& _positions;
        kmer_index<alphabet_t, k> _index;
};