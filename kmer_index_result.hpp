// Copyright (c) 2020 Clemens Cords. All rights reserved.

#pragma once

#include <kmer_index.hpp>
#include <compressed_bitset.hpp>

namespace kmer::detail
{
    // forward declaration
    template<seqan3::alphabet, typename, size_t>
    class kmer_index_element;

    // result type that only holds pointers to the positions inside kmer index
    // only to be used internally
    template<typename position_t>
    struct kmer_index_result
    {
        private:
            // iterator class that automatically jumps to next valid result
            class kmer_index_result_iterator
            {
                friend class kmer_index_result<position_t>;

                private:
                    // result the iterator is operating on
                    const kmer_index_result<position_t>* _result;

                    // current index of the position pointed to
                    size_t _position_i = 0;

                    // first and last valid index for the positions, computer on construction
                    size_t _first_valid_i, _last_valid_i;

                    // move iterator to next valid position
                    // returns true if moving possible, false if end is hit
                    bool advance_to_next_valid_result()
                    {
                        if (_position_i == _last_valid_i)
                            return false;

                        size_t new_i = _position_i + 1;

                        while (new_i < _result->_bitmask.size() and _result->_bitmask.at(new_i) == false)
                            new_i++;

                        assert(new_i != _position_i);

                        _position_i = new_i;
                        return true;
                    }

                    // move iterator to previous valid position
                    // returns true if moving possible, false if start is hit
                    bool advance_to_previous_valid_result()
                    {
                        if (_position_i == _first_valid_i)
                            return false;

                        size_t new_i = _position_i - 1;

                        while (new_i > 0 and _result->bitmask.at(new_i) == false)
                            new_i--;

                        assert(new_i != _position_i);

                        _position_i = new_i;
                        return true;
                    }

                protected:
                    // ctor, choose wether to start at beginning or end of range
                    kmer_index_result_iterator(kmer_index_result<position_t>* result, bool beginning_or_end)
                            : kmer_index_result_iterator(result)
                    {
                        if (beginning_or_end)
                        {
                            if (not _result->_bitmask.at(0))
                                advance_to_next_valid_result();
                        }
                        else
                        {
                            _position_i = _last_valid_i+1;
                        }
                    }

                public:
                    // typedefs required by std
                    using iterator_category = std::bidirectional_iterator_tag;
                    using value_type = position_t;
                    using difference_type = void;
                    using pointer = void;
                    using reference = void;

                    // typedef of own type for readability
                    using iterator_t = kmer_index_result<position_t>::kmer_index_result_iterator;

                    // CTOR
                    kmer_index_result_iterator(kmer_index_result<position_t>* result)
                            : _result(result), _position_i(0)
                    {
                        size_t i = result->_bitmask.size() - 1;
                        for (i; i >= 0; i--)
                            if (result->_bitmask.at(i))
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

                    // operators: arithmetics
                    iterator_t& operator++()
                    {
                        advance_to_next_valid_result();
                        return *this;
                    }

                    iterator_t& operator+=(int i)
                    {
                        if (i > 0)
                        {
                            while (i > 0)
                            {
                                advance_to_next_valid_result();
                                i--;
                            }
                            return *this;
                        }
                        else
                        {
                            while (i < 0)
                            {
                                advance_to_previous_valid_result();
                                i++;
                            }
                            return *this;
                        }
                    }

                    iterator_t& operator--()
                    {
                        advance_to_previous_valid_result();
                        return *this;
                    }

                    iterator_t& operator-=(int i)
                    {
                        if (i > 0)
                        {
                            while (i > 0)
                            {
                                advance_to_previous_valid_result();
                                i--;
                            }
                            return *this;
                        }
                        else
                        {
                            while (i < 0)
                            {
                                advance_to_next_valid_result();
                                i++;
                            }
                            return this;
                        }
                    }

                    // compare
                    bool operator==(iterator_t other)
                    {
                        return this->_position_i == other._position_i and this->_result == other._result;
                    }

                    bool operator!=(iterator_t other)
                    {
                        return not(*this == other);
                    }

                    // dereference
                    value_type operator*()
                    {
                        auto result = _result->at(_position_i);
                        return result;
                    }
            };

            // bitmask specifies which of the results should be ignore
            compressed_bitset<uint_fast64_t> _bitmask;

            // pointers to positions inside kmer index map
            std::vector<const std::vector<position_t>*> _positions;

        public:
            // at, throws if the position at i isn't usable
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

            // CTOR: default
            explicit kmer_index_result(bool zero_or_one = true)
                    : _bitmask(0, zero_or_one)
            {
            }

            // CTOR: used by nk and k search
            explicit kmer_index_result(const std::vector<position_t>* positions, bool zero_or_one = true)
                    : _bitmask(positions->size(), zero_or_one), _positions{positions}
            {
            }

            // CTOR: used by subk searchs
            explicit kmer_index_result(std::vector<const std::vector<position_t>*> positions, bool zero_or_one = true)
                    : _positions(positions.begin(), positions.end()),
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

            // get number of valid positions
            size_t size() const
            {
                return _bitmask.count_bits_equal_to(true);
            }

            // copy result into vector and sort
            // terrible runtime
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

                if (sort_results)
                    std::sort(output.begin(), output.end());

                return output;
            }

            // iterables
            kmer_index_result_iterator begin()
            {
                return kmer_index_result_iterator(this, true); // set to beginning
            }

            kmer_index_result_iterator end()
            {
                return kmer_index_result_iterator(this, false); // set to end
            }
    };

} // end of namespace kmer::detail
