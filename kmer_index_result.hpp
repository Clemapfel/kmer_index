// Copyright (c) 2020 Clemens Cords. All rights reserved.

#pragma once

#include <kmer_index.hpp>
#include <bitset>

namespace detail {
    // runtime-optimized equivalent for std::vector<bool>
    // is created in (size / sizeof(integer_t)) rather than size like std::vector<bool>(size)
    template<typename integer_t>
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
                    : _bits(n_bits & _and_v, _not_zero) {
            }

            // set ith bit to 0
            void set_0(size_t i) {
                assert(i < n_bits);

                size_t n = i & _and_v;
                _bits.data()[i >> _rshift_v] |= _zero << n;
            }

            // set ith bit to 1
            void set_1(size_t i) {
                assert(i < n_bits);

                size_t n = i & _and_v;
                _bits.data()[i >> _rshift_v] |= _one << n;
            }

            // get ith bit
            bool at(size_t i) {
                assert(i < n_bits);

                size_t n = i & _and_v;
                return (_bits[i >> _rshift_v] & (_one << n)) != _zero;
            }

            // resets all bits to 1
            void clear() {
                for (auto& i : _bits)
                    i = ~_zero;
            }
    };

// result proxy that supports lazy-eval get
    template<seqan3::alphabet alphabet_t, size_t k, typename position_t = uint32_t>
    class kmer_index_results {
            friend class kmer_index<alphabet_t, k>;

        public:

        protected:
            explicit kmer_index_results(const std::vector<position_t>& positions)
                    : _bitmask(positions.size()) {
                _positions = positions;
                //_bitmask.reserve(positions.size());
            }

            void set_should_use(size_t i) {
                _bitmask.set_0(i);
            }

        private:
            const std::vector<position_t>& _positions;
            compressed_bitset<size_t> _bitmask;
    };
}