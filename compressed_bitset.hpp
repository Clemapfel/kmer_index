// Copyright (c) 2020 Clemens Cords. All rights reserved.

#pragma once


namespace kmer
{
namespace detail
{
    // runtime-optimized equivalent for std::vector<bool>
    template<typename integer_t = uint_fast64_t>
    class compressed_bitset
    {
        private:
            // pre-calculate frequently used constants
            static constexpr integer_t _and_v = (sizeof(integer_t) * 8) - 1;          // i % n = i & n-1
            static constexpr integer_t _rshift_v = log2((sizeof(integer_t) * 8));  // i / n = i >> log2(n)
            static constexpr integer_t _one = 1, _zero = 0, _not_zero = ~_zero;

            // last valid index (equivalent to std::vector::size())
            size_t _n_bits;

            // holds integers
            std::vector<integer_t> _bits;

        public:
            // CTOR
            // while n_bits is arbitrary, if n_bits % sizeof(integer_t) != 0, there will be overfill
            compressed_bitset(size_t n_bits, bool zero_or_one)
                    : _bits(std::max(n_bits / (sizeof(integer_t) * 8) + 1, 1ul), (zero_or_one ? _not_zero : _zero))
            {
                _n_bits = n_bits;
            }

            // convert to regular vector
            std::vector<bool> to_vector() const
            {
                std::vector<bool> out;
                for (size_t i = 0; i < _n_bits; ++i)
                    out.push_back(at(i));

                return out;
            }

            // cast to regular vector
            explicit operator std::vector<bool>()
            {
                return to_vector();
            }


            // set ith bit to 0
            void set_0(size_t i)
            {
                if (i >= _n_bits)
                    throw std::out_of_range("compressed bitset index out of range");

                size_t n = i & _and_v;
                _bits.data()[i >> _rshift_v] &= ~(_one << n);
            }

            // set ith bit to 1
            void set_1(size_t i)
            {
                if (i >= _n_bits)
                    throw std::out_of_range("compressed bitset index out of range");

                size_t n = i & _and_v;
                _bits.data()[i >> _rshift_v] |= _one << n;
            }

            // get ith bit
            bool at(size_t i) const
            {
                if (i >= _n_bits)
                    throw std::out_of_range("compressed bitset index out of range");

                size_t n = i & _and_v;
                auto v = (_bits[i >> _rshift_v] & (_one << n)) != _zero;
                return v;
            }

            // reset all bits to 1
            void clear_to_1()
            {
                for (auto& i : _bits)
                    i = _not_zero;
            }

            // reset all bits to 0
            void clear_to_0()
            {
                for (auto& i : _bits)
                    i = _zero;
            }

            size_t size() const
            {
                return _n_bits;
            }

            // popcount
            size_t count_bits_equal_to(bool b) const
            {
                // TODO: do new x86 std::popcount
                size_t n_ones = 0;

                for (size_t i = 0; i < _n_bits; ++i)
                    if (at(i))
                        n_ones++;

                return (b ? n_ones : _n_bits - n_ones);
            }
    };

} // end of namespace detail
} // end of namespace kmer