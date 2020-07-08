// Copyright (c) 2020 Clemens Cords. All rights reserved.

#pragma once

#include <kmer_index.hpp>
#include <compressed_bitset.hpp>

namespace kmer::detail
{
    enum class BYPASS_BITMASK : bool {YES = true, NO = false};

    template<typename position_t>
    class kmer_index_result
    {
        private:
            compressed_bitset<uint_fast64_t> _bitmask;
            const bool _bypass_bitmask;
            const size_t _n_results;

            // pointers to positions inside kmer index map
            std::vector<const std::vector<position_t> *> _positions;

        protected:
            // iterator class that automatically jumps to next valid result
            class kmer_index_result_iterator
            {
                friend class kmer_index_result<position_t>;

                private:
                    // result the iterator is operating on
                    const kmer_index_result<position_t> *_result;

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

                        while (new_i < _result->get_n_results() and _result->is_valid(new_i) == false)
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

                        while (new_i > 0 and _result->is_valid(new_i) == false)
                            new_i--;

                        assert(new_i != _position_i);

                        _position_i = new_i;
                        return true;
                    }

                protected:
                    // ctor, choose wether to start at beginning or end of range
                    kmer_index_result_iterator(kmer_index_result<position_t> *result, bool beginning_or_end)
                            : kmer_index_result_iterator(result)
                    {
                        if (beginning_or_end)
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
                    // typedefs required by std
                    using iterator_category = std::bidirectional_iterator_tag;
                    using value_type = position_t;
                    using difference_type = void;
                    using pointer = void;
                    using reference = void;

                    // typedef of own type for readability
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

                    // operators: arithmetics
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
            // CTOR: defailt
            kmer_index_result()
                    : _bitmask(0, true), _bypass_bitmask(false), _n_results(0)
            {
            }

            // CTOR: k or nk
            kmer_index_result(const std::vector<position_t> *positions, bool zero_or_one, BYPASS_BITMASK bypass_bitmask)
                    : _bitmask((bool(bypass_bitmask) ? 0 : positions->size()), zero_or_one), _bypass_bitmask(bool(bypass_bitmask)), _n_results(positions->size())
            {
                _positions = {positions};
            }

            // CTOR: sub k
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

            // get number of valid positions
            size_t size() const
            {
                return _bitmask.count_bits_equal_to(true);
            }

            // copy result into vector and sort
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



