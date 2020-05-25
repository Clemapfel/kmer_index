// Copyright (c) 2020 Clemens Cords. All rights reserved.

#pragma once

#include <kmer_index.hpp>

template<seqan3::alphabet, size_t, typename>
class kmer_index_element;

namespace detail
{
    // runtime-optimized equivalent for std::vector<bool>
    template<typename integer_t = uint_fast64_t>
    class compressed_bitset{
        private:
            // pre-calculate frequently used constants
            static constexpr integer_t _and_v = (sizeof(integer_t) * 8) - 1;          // i % n = i & n-1
            static constexpr integer_t _rshift_v = log2((sizeof(integer_t) * 8));  // i / n = i >> log2(n)
            static constexpr integer_t _one = 1, _zero = 0, _not_zero = ~_zero;

            size_t _n_bits;

            // holds integers
            std::vector<integer_t> _bits;

        public:
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

            size_t size() const
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
    struct kmer_index_result
    {
       friend class kmer_index_element<alphabet_t, k, position_t>;

        private:
            // iterator class that automatically jumps to next valid result
            class kmer_index_result_iterator
            {
                friend class kmer_index_result<alphabet_t, k, position_t>;

                private:
                    const kmer_index_result<alphabet_t, k, position_t>* _result;
                    size_t _position_i = 0;
                    size_t _first_valid_i, _last_valid_i;

                    bool advance_to_next_valid_result()
                    {
                        if (_position_i == _last_valid_i)
                            return false;

                        size_t new_i = _position_i+1;

                        while(new_i < _result->_bitmask.size() and _result->_bitmask.at(new_i) == false)
                            new_i++;

                        assert(new_i != _position_i);

                        _position_i = new_i;
                        return true;
                    }

                    bool advance_to_previous_valid_result()
                    {
                        if (_position_i == _first_valid_i)
                            return false;

                        size_t new_i = _position_i-1;

                        while (new_i > 0 and _result->bitmask.at(new_i) == false)
                            new_i--;

                        assert(new_i != _position_i);

                        _position_i = new_i;
                        return true;
                    }

                protected:
                    kmer_index_result_iterator(kmer_index_result<alphabet_t, k, position_t>* result, bool beginning_or_end)
                        : kmer_index_result_iterator(result)
                    {
                        if (beginning_or_end)
                        {
                            if (not _result->_bitmask.at(0))
                                advance_to_next_valid_result();
                        }
                        else
                        {
                            _position_i = _last_valid_i;
                        }
                    }

                public:
                    // typedefs required
                    using iterator_category = std::bidirectional_iterator_tag;
                    using value_type = position_t;
                    using difference_type = void;
                    using pointer = void;
                    using reference = void;

                    using iterator_t = kmer_index_result<alphabet_t, k, position_t>::kmer_index_result_iterator;

                    // ctor
                    kmer_index_result_iterator(kmer_index_result<alphabet_t, k, position_t>* result)
                        : _result(result), _position_i(0)
                    {
                        size_t i= result->_bitmask.size()-1;
                        for (i; i >= 0; i--)
                            if(result->_bitmask.at(i))
                                break;

                        _last_valid_i = i;

                        if (not result->_bitmask.at(0))
                        {
                            advance_to_next_valid_result();
                            _first_valid_i = _position_i;
                        }
                        else
                            _first_valid_i = 0;
                    }

                    // operators
                    iterator_t& operator++()
                    {
                        advance_to_next_valid_result();
                        return *this;
                    }

                    iterator_t& operator++(int i)
                    {
                        while(i > 0)
                        {
                            advance_to_next_valid_result();
                            if (_position_i == _result->size())
                                break;

                            i--;
                        }
                        return *this;
                    }

                    iterator_t& operator--()
                    {
                        advance_to_previous_valid_result();
                        return *this;
                    }

                    iterator_t& operator--(int i)
                    {
                        while (i > 0)
                        {
                            advance_to_previous_valid_result();
                            if (_position_i == 0)
                                break;

                            i--;
                        }
                        return *this;
                    }

                    bool operator==(iterator_t other)
                    {
                        return this->_position_i == other._position_i and this->_result == other._result;
                    }

                    bool operator!=(iterator_t other)
                    {
                        return not (*this == other);
                    }

                    value_type operator*()
                    {
                        auto result = _result->at(_position_i);
                        return result;
                    }
            };

            using index_t = kmer_index_element<alphabet_t, k, position_t>;

            // bitmask specifies which of the results should be ignore
            compressed_bitset<uint_fast64_t> _bitmask;

            // pointers to positions inside kmer index map
            const std::vector<const std::vector<position_t>*> _positions;

            //keep index in memory so results don't become invalid
            const index_t& _index;

        protected:
            position_t at(size_t i) const
            {
                assert(_bitmask.at(i));

                if (i < _positions.at(0)->size())
                    return _positions.at(0)->at(i);

                size_t vector_i = 0;
                for (; vector_i < _positions.size(); ++vector_i)
                {
                    if (int(i) - int(_positions.at(vector_i)->size()) < 0)
                        break;

                    else
                        i -= _positions.at(vector_i)->size();
                }

                return _positions.at(vector_i)->at(i);
            }

            // ctors only used by kmer_index
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

            const std::vector<const std::vector<position_t>*>* get_positions_raw() const
            {
                return &_positions;
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
            // user should not be able to construct results, only kmer_index does
            //kmer_index_result() = default;
            //kmer_index_result(const kmer_index_result&) = delete;
            //kmer_index_result(kmer_index_result&&) = delete;

            // get number of valid positions
            size_t size() const
            {
                return _bitmask.count_bits_equal_to(true);
            }

            size_t count() const
            {
                return size();
            }

            // copy result into vector
            std::vector<position_t> to_vector(bool sort_results = true) const
            {
                if (_positions.empty())
                    return std::vector<position_t>();

                std::vector<position_t> output;

                size_t i = 0;
                for (const auto* vec : _positions)
                    for (size_t j = 0; j < vec->size(); ++j, ++i)
                        if (_bitmask.at(i))
                            output.push_back(vec->at(j));

                if(sort_results)
                    std::sort(output.begin(), output.end());

                return output;
            }

            // iterables
            kmer_index_result_iterator begin()
            {
                return kmer_index_result_iterator(this, true);
            }

            kmer_index_result_iterator end()
            {
                return kmer_index_result_iterator(this, false);
            }
    };
} // end of namespace detail