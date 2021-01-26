// Copyright (c) 2020 Clemens Cords. All rights reserved.

#pragma once

#include <kmer_index.hpp>
#include <compressed_bitset.hpp>

namespace kmer::detail
{
    // if bitmask bypassed, skip operator arithmetics and treat bitmask as all 11111...11
    enum class BYPASS_BITMASK : bool {YES = true, NO = false};

    // container for kmer index results (c.f. [1])
    template<typename position_t>
    class kmer_index_result
    {
        private:
            compressed_bitset<uint_fast64_t> _bitmask;
            const bool _bypass_bitmask;
            const size_t _n_results;

            // pointers to positions inside kmer_index::_data
            std::vector<const std::vector<position_t>*> _positions;

        protected:
            class kmer_index_result_iterator
            {
                friend class kmer_index_result<position_t>;

                private:
                    const kmer_index_result<position_t> *_result;
                    size_t _position_i = 0;

                    // first and last valid index for the positions, computed at construction
                    size_t _first_valid_i, _last_valid_i;

                    bool advance_to_next_valid_result()
                    {
                        if (_position_i == _last_valid_i)
                            return false;

                        size_t new_i = _position_i + 1;

                        while (new_i < _result->get_n_results() and _result->is_valid(new_i) == false)
                            new_i++;

                        assert(new_i != _position_i);

                        _position_i = new_i;
                        return true;
                    }

                    bool advance_to_previous_valid_result()
                    {
                        if (_position_i == _first_valid_i)
                            return false;

                        size_t new_i = _position_i - 1;

                        while (new_i > 0 and _result->is_valid(new_i) == false)
                            new_i--;

                        assert(new_i != _position_i);

                        _position_i = new_i;
                        return true;
                    }

                protected:
                    kmer_index_result_iterator(kmer_index_result<position_t> *result, bool start_at_beginning_or_end)
                            : kmer_index_result_iterator(result)
                    {
                        if (start_at_beginning_or_end)
                        {
                            if (not _result->_bitmask.at(0))
                                advance_to_next_valid_result();
                        }
                        else
                        {
                            _position_i = _last_valid_i + 1;
                        }
                    }

                public:
                    using iterator_category = std::bidirectional_iterator_tag;
                    using value_type = position_t;
                    using difference_type = void;
                    using pointer = void;
                    using reference = void;

                    // typedef for readability
                    using iterator_t = kmer_index_result<position_t>::kmer_index_result_iterator;

                    // CTOR
                    kmer_index_result_iterator(kmer_index_result<position_t> *result)
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

                    iterator_t &operator++()
                    {
                        advance_to_next_valid_result();
                        return *this;
                    }

                    iterator_t &operator+=(int i)
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

                    iterator_t &operator--()
                    {
                        advance_to_previous_valid_result();
                        return *this;
                    }

                    iterator_t &operator-=(int i)
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

                    bool operator==(iterator_t other)
                    {
                        return this->_position_i == other._position_i and this->_result == other._result;
                    }

                    bool operator!=(iterator_t other)
                    {
                        return not(*this == other);
                    }

                    value_type operator*()
                    {
                        auto result = _result->at(_position_i);
                        return result;
                    }
            };

            // is i a valid position
            bool is_valid(size_t i) const
            {
                if (_bypass_bitmask)
                    return true;
                else
                    return _bitmask.at(i);
            }

            size_t get_n_results() const
            {
                return _n_results;
            }

        public:
            // CTORs
            kmer_index_result()
                    : _bitmask(0, true), _bypass_bitmask(false), _n_results(0)
            {
            }

            kmer_index_result(const std::vector<position_t> *positions, bool fill_with_zero_or_ones, BYPASS_BITMASK bypass_bitmask)
                    : _bitmask((bool(bypass_bitmask) ? 0 : positions->size()), fill_with_zero_or_ones), _bypass_bitmask(bool(bypass_bitmask)), _n_results(positions->size())
            {
                _positions = {positions};
            }

            kmer_index_result(std::vector<const std::vector<position_t>*> positions)
                    : _bitmask(0, true),
                      _bypass_bitmask(true),
                      _n_results([&]() -> size_t {
                            size_t n = 0;
                            for (const auto vec : positions)
                                n += vec->size();
                            return n;
                        }())
            {
                _positions = positions;
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

            // number of valid positions
            size_t size() const
            {
                return _bitmask.count_bits_equal_to(true);
            }

            std::vector<position_t> to_vector() const
            {
                if (_positions.empty())
                    return std::vector<position_t>();

                std::vector<position_t> output;

                size_t i = 0;
                for (const auto *vec : _positions)
                    for (size_t j = 0; j < vec->size(); ++j, ++i)
                        if (is_valid(i))
                            output.push_back(vec->at(j));

                std::sort(output.begin(), output.end());

                return output;
            }

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

// ###################################
//
// [1]
//
// When searching for query positions with a query of length != k, the kmer index has to not
// only lookup the corresponding vector of positions in it's data member but furthermore check
// each position for wether the query actually occurs there. Usually one would just discard
// the non-valid positions however this requires reallocation. This kmer-index result class allows
// the positional vector (holding all positions of the k-prefix of the quey) to be returned by
// reference by also holding a bitmask that for each index in that vector marks wether or not the
// rest of the query actually occurs after it.
// By utilizing kmer_index_result instead of simple std::vector and furthermore optimizing the bitmask
// to be compressed the minimal amount of allocation per query is needed drastically improving performance
//
// ###################################



