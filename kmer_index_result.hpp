// Copyright (c) 2020 Clemens Cords. All rights reserved.

#pragma once

#include <kmer_index.hpp>
#include <bitset>

template<seqan3::alphabet alphabet_t, size_t k, typename position_t = uint32_t>
class kmer_index_results
{
    friend class kmer_index<alphabet_t, k>;

    public:

    protected:
        kmer_index_results(const std::vector<position_t>& positions)
        {
            _positions = positions;
            _bitmask.reserve(_positions.size());

            for (size_t i = 0; i < _positions.size(); ++i)
                _bitmask.emplace_back(true);
        }

        void set_should_use(size_t i)
        {
            _bitmask.data()[i] = false;
        }

    private:
        const std::vector<position_t>& _positions;
        std::vector<bool> _bitmask;
        kmer_index<alphabet_t, k> _index;
};