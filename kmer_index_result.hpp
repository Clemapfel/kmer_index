// Copyright (c) 2020 Clemens Cords. All rights reserved.

#pragma once

#include <kmer_index.hpp>

template<seqan3::alphabet, size_t, typename>
class kmer_index_element;

namespace detail
{
    // runtime-optimized equivalent for std::vector<bool>
    template<typename integer_t = uint_fast64_t>
    class compressed_bitset {
        private:
            // pre-calculate frequently used constants
            static constexpr integer_t _and_v = (sizeof(integer_t) * 8) - 1;          // i % n = i & n-1
            static constexpr integer_t _rshift_v = log2((sizeof(integer_t) * 8));  // i / n = i >> log2(n)
            static constexpr integer_t _one = 1, _zero = 0, _not_zero = ~_zero;

            size_t _n_bits;

            // holds integers
            std::vector<integer_t> _bits;

        protected:
            // convert to regular vector
            std::vector<bool> to_vector() const {
                std::vector<bool> out;
                for (size_t i = 0; i < _n_bits; ++i)
                    out.push_back(at(i));

                return out;
            }

        public:
            // create by specifying maximum number of bits
            // while supporting an arbitrary number, unless n_bits & sizeof(integer_t) == 0 more bits than necessary
            // have to be allocated. Consider using a smaller integer_t if this proofs problematic
            compressed_bitset(size_t n_bits, bool zero_or_one)
                    : _bits(std::max(n_bits / (sizeof(integer_t) * 8), 1ul), (zero_or_one ? _not_zero : _zero))
            {
                _n_bits = n_bits;
            }

            // set ith bit to 0
            void set_0(size_t i)
            {
                assert(i < _n_bits);

                size_t n = i & _and_v;
                _bits.data()[i >> _rshift_v] &= ~(_one << n);
            }

            // set ith bit to 1
            void set_1(size_t i)
            {
                assert(i < _n_bits);

                size_t n = i & _and_v;
                _bits.data()[i >> _rshift_v] |= _one << n;
            }

            // get ith bit
            bool at(size_t i) const
            {
                assert(i < _n_bits);

                size_t n = i & _and_v;
                return (_bits[i >> _rshift_v] & (_one << n)) != _zero;
            }

            // resets all bits to 1
            void clear_to_1()
            {
                for (auto& i : _bits)
                    i = _not_zero;
            }

            // to 0
            void clear_to_0()
            {
                for (auto& i : _bits)
                    i = _zero;
            }

            size_t size()
            {
                return _n_bits;
            }

            size_t count_bits_equal_to(bool b) const
            {
                size_t n_ones = 0;

                for (size_t i = 0; i < _n_bits; ++i)
                    if (at(i))
                        n_ones++;

                return (b ? n_ones : _n_bits - n_ones);
            }

            // cast to regular vector
            explicit operator std::vector<bool>() {
                return to_vector();
            }
    };

    // result type that only holds pointers to the positions inside kmer index
    template<seqan3::alphabet alphabet_t, size_t k, typename position_t>
    class kmer_index_result
    {
       friend class kmer_index_element<alphabet_t, k, position_t>;

        private:
            using index_t = kmer_index_element<alphabet_t, k, position_t>;

            // bitmask specifies which of the results should be ignore
            compressed_bitset<uint_fast64_t> _bitmask;

            // pointers to positions inside kmer index map
            const std::vector<const std::vector<position_t>*> _positions;

            //keep index in memory so results don't become invalid
            const index_t& _index;

        //protected:
        public:
            // ctors
            kmer_index_result(const index_t* index, bool zero_or_one = true)
                    : _index(*index), _bitmask(0, zero_or_one)
            {
            }

            kmer_index_result(const std::vector<position_t>* positions, const index_t* index, bool zero_or_one = false)
                    : _index(*index), _bitmask(positions->size(), zero_or_one), _positions{positions}
            {

            }

            kmer_index_result(std::vector<const std::vector<position_t>*> positions, const index_t* index, bool zero_or_one = false)
                    : _index(*index),
                      _positions(positions.begin(), positions.end()),
                      _bitmask([&positions]() {
                          size_t n = 0;
                          for (const auto* p : positions)
                              n += p->size();
                          return n;
                      }(), zero_or_one)
            {
            }

            // specify which positions to use by setting bitmask
            void should_not_use(size_t i)
            {
                _bitmask.set_0(i);
            }

            void should_use(size_t i)
            {
                _bitmask.set_1(i);
            }

        public:
            size_t size() const
            {
                return _bitmask.count_bits_equal_to(true);
            }

            // lazy eval placeholder
            std::vector<position_t> to_vector() const
            {
                if (_positions.empty())
                    return std::vector<position_t>();

                std::vector<position_t> output;

                size_t i = 0;
                for (const auto* vec : _positions)
                    for (size_t j = 0; j < vec->size(); ++j, ++i)
                        if (_bitmask.at(i))
                            output.push_back(vec->at(j));

                std::sort(output.begin(), output.end());
                return output;
            }

    };
}