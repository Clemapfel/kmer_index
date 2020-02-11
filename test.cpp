//
// Created by Clemens Cords on 2/7/20.
// Copyright (c) 2020 Clemens Cords. All rights reserved.
//

#include <gtest/gtest.h>

#include <kmer_index.hpp>
#include <input_generator.hpp>

#include <seqan3/search/fm_index/fm_index.hpp>
#include <seqan3/alphabet/concept.hpp>

template<size_t text_length, size_t n_queries>
auto compare_with_fm()
{
    auto text = input_generator<dna5, 






}
