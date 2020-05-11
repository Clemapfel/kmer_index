// Copyright (c) 2020 Clemens Cords. All rights reserved.

#pragma once

#include <kmer_index.hpp>
#include <bitset>

namespace detail {
    // runtime-optimized equivalent for std::vector<bool>
    // is created in (size / sizeof(integer_t)) rather than size like std::vector<bool>(size)
    template<typename integer_t = uint_fast64_t>
    class compressed_bitset {
        private:
            // pre-calculate frequently used constants
            static constexpr integer_t _and_v = sizeof(integer_t) - 1;       // i % n = i & n-1
            static constexpr integer_t _rshift_v = log2(sizeof(integer_t)); // i / n = i >> log2(n)

            static constexpr integer_t _one = 1, _zero = 0, _not_zero = ~_zero;

            // holds integers
            std::vector<integer_t> _bits;

        public:
            // create by specifying maximum i
            compressed_bitset(size_t n_bits)
                    : _bits(n_bits & _and_v, _not_zero)
            {
            }

            // set ith bit to 0
            void set_0(size_t i)
            {
                assert(i < n_bits);

                size_t n = i & _and_v;
                _bits.data()[i >> _rshift_v] |= _zero << n;
            }

            // set ith bit to 1
            void set_1(size_t i)
            {
                assert(i < n_bits);

                size_t n = i & _and_v;
                _bits.data()[i >> _rshift_v] |= _one << n;
            }

            // get ith bit
            bool at(size_t i) const
            {
                assert(i < n_bits);

                size_t n = i & _and_v;
                return (_bits[i >> _rshift_v] & (_one << n)) != _zero;
            }

            // resets all bits to 1
            void clear()
            {
                for (auto& i : _bits)
                    i = ~_zero;
            }
    };
}

// result proxy that supports lazy-eval get
    template<seqan3::alphabet alphabet_t, size_t k, typename position_t = uint32_t>
    class kmer_index_results
    {
        friend class kmer_index<alphabet_t, k>;

        public:
            std::vector<position_t> to_vector() const
            {
                std::vector<position_t> output;

                size_t i = 0;
                for (const auto* vec : _positions)
                    for (size_t j = 0; j < vec->size(); ++j, ++i)
                        if (_bitmask.at(i))
                            output.push_back(vec->at(j));

                return output;
            }

        //protected:
            explicit kmer_index_results(const std::vector<position_t>* positions)
                    : _bitmask(positions->size())
            {
                _positions.reserve(1);
                _positions.push_back(positions);
            }

            explicit kmer_index_results(std::vector<const std::vector<position_t>*> positions)
                    : _positions(positions.begin(), positions.end()),
                      _bitmask([positions]() {
                          size_t n = 0;
                          for (const auto* p : positions)
                            n += p->size();
                          return n;
                      }())
            {

            }

            void set_should_use(size_t i) {
                _bitmask.set_0(i);
            }

        private:
            std::vector<const std::vector<position_t>*> _positions;
            detail::compressed_bitset<size_t> _bitmask;
    };