// ============================================================================
//                 SeqAn - The Library for Sequence Analysis
// ============================================================================
//
// Copyright (c) 2006-2017, Knut Reinert & Freie Universitaet Berlin
// Copyright (c) 2016-2017, Knut Reinert & MPI Molekulare Genetik
// All rights reserved.
//
// Redistribution and use in source and binary forms, with or without
// modification, are permitted provided that the following conditions are met:
//
//     * Redistributions of source code must retain the above copyright
//       notice, this list of conditions and the following disclaimer.
//     * Redistributions in binary form must reproduce the above copyright
//       notice, this list of conditions and the following disclaimer in the
//       documentation and/or other materials provided with the distribution.
//     * Neither the name of Knut Reinert or the FU Berlin nor the names of
//       its contributors may be used to endorse or promote products derived
//       from this software without specific prior written permission.
//
// THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
// AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
// IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
// ARE DISCLAIMED. IN NO EVENT SHALL KNUT REINERT OR THE FU BERLIN BE LIABLE
// FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL
// DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR
// SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER
// CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT
// LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY
// OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH
// DAMAGE.
//
// ============================================================================

//! \cond DEV

/*! \file
 *  \brief Provides seqn3::aligned_sequence_adaptor_constant_access.
 *  \ingroup container
 *  \author Marie Hoffmann <marie.hoffmann AT fu-berlin.de>
 */

#pragma once

#include <algorithm>
#include <initializer_list>
#include <iostream>
#include <string>
#include <type_traits>
#include <vector>

#include <seqan3/alphabet/concept.hpp>
#include <seqan3/alphabet/gap/gapped.hpp>
#include <seqan3/range/container/concept.hpp>
#include <seqan3/range/detail/random_access_iterator.hpp>

#include "sdsl/bit_vectors.hpp"
#include "sdsl/rank_support_v5.hpp"
#include "sdsl/select_support_mcl.hpp"
#include "sdsl/util.hpp"

namespace seqan3 {

/*!\brief Container for an aligned sequence with random access at constant time.
 * \tparam inner_type {The container type with which aligned sequence structure
 * can be constructed, e.g. std::vector<gapped<dna4>>. The inner type of the
 * container must fulfill the seqan3::alphabet_concept.}
 * \implements partially the seqan3::random_access_sequence_concept
 * \remark {This class has to be parameterized with the union type
 * gapped<alphabet_type> packed into a container supporting the random access
 * operators [] and at(). This way an aligned sequence can be initialized or
 * assigned to any sequence of alphabet and/or gap symbols.}
 *
 * An aligned sequence is internally represented by a bit_vector
 * set to one whenever an aligned position corresponds to a gap. The ungapped
 * sequence making up the 0s in the bit vector is stored sequentially in a vector.
 * We use a rank and a select function on the bit vector to map from the gapped
 * sequence space to the ungapped sequence space and vice versa.
 * \par Example
 * '---AT--ATC-GT' is stored as a bit vector '1110011000100' and a compressed
 * sequence 'ATATCGT'.
 */

 template <typename type>
 concept bool random_access_concept = requires (type val)
 {
     // member types
     typename type::reference;
     // access container
     { val[0]    } -> typename type::reference;
     { val.at(0) } -> typename type::reference;
 };

template <typename inner_type>
//!\cond
    //requires random_access_sequence_concept<inner_t> && alphabet_concept<ranges::v3::value_type_t<inner_type>>
    requires random_access_concept<inner_type> &&
    alphabet_concept<ranges::v3::value_type_t<inner_type>>
//!\endcond
struct aligned_sequence_adaptor_constant_access
{

private:
     //!\privatesection
    using aligned_sequence_t = aligned_sequence_adaptor_constant_access;
    //using alphabet_t = typename ranges::v3::value_type_t<inner_type>::alphabet_type;

// re-initialize rank support structure of bit_vector
void update_support_structures()
{
    rs = sdsl::rank_support_v5<1, 1>(&gap_vector);
    ss = sdsl::select_support_mcl<0, 1>(&gap_vector);
}

public:
    //!\publicsection
    /*!\name Member types
    * \{
    */
    //!\brief Value type of container elements.
    //!\hideinitializer
    using value_type = typename ranges::v3::value_type_t<inner_type>;

    //!\brief Use reference type defined by container.
    //!\hideinitializer
    using reference = value_type;

    //!\brief Use const reference type provided by container.
    //!\hideinitializer
    using const_reference = const reference;

    //!\brief Use random access iterator on container as iterator type.
    //!\hideinitializer
    using iterator = detail::random_access_iterator<aligned_sequence_t>;

    //!\brief Use const random access iterator on container as const iterator type.
    //!\hideinitializer
    using const_iterator = detail::random_access_iterator<aligned_sequence_t const>;

    //!\brief Type for distances between iterators is taken from alphabet container.
    //!\hideinitializer
    using difference_type = typename ranges::v3::difference_type_t<inner_type>;

    //!\brief Use alphabet container's size_type as a position.
    //!\hideinitializer
    using size_type = typename ranges::v3::size_type_t<inner_type>;
    //!\}

    /* rule of six */
    /*!\name Constructors, destructor and assignment
    * \{
    */
    // \brief Default constructor.
    constexpr aligned_sequence_adaptor_constant_access() = default;

    //!\brief Default copy constructor.
    constexpr aligned_sequence_adaptor_constant_access(aligned_sequence_adaptor_constant_access const &) = default;

    //!\brief Default copy construction via assignment.
    constexpr aligned_sequence_adaptor_constant_access & operator=(aligned_sequence_adaptor_constant_access const &) = default;

    /*!\brief Move constructor.
    *
    * The rank support structure of the bit vector needs the moved reference to
    * the latter one. It has to be hand over subsequently. Therefore no default
    * move assignment is possible. See https://github.com/xxsds/sdsl-lite/issues/15.
    */
    constexpr aligned_sequence_adaptor_constant_access (aligned_sequence_adaptor_constant_access && rhs) :
    gap_vector(std::move(rhs.gap_vector)), sequence(std::move(rhs.sequence))
    {
        rs.set_vector(&gap_vector);
        ss.set_vector(&gap_vector);
    }

    //!\brief Move assignment.
    constexpr aligned_sequence_adaptor_constant_access & operator=(aligned_sequence_adaptor_constant_access && rhs)
    {
        gap_vector = std::move(rhs.gap_vector);
        sequence = std::move(rhs.sequence);
        rs.set_vector(&gap_vector);
        ss.set_vector(&gap_vector);
        return *this;
    }

    //!\brief Use default deconstructor.
    ~aligned_sequence_adaptor_constant_access() = default;
    //!\}

    //!\brief
    /*!\name Constructors of sequence concept.
    * \{
    */
    //!\brief Construct by single value repeated 'size' times.
    constexpr aligned_sequence_adaptor_constant_access(size_type size, value_type value)
    {
        assign(size, value);
    };

    //!\brief Construct by range given by two iterators of an aligned sequence.
    constexpr aligned_sequence_adaptor_constant_access(iterator it1, iterator it2)
    {
        assign(it1, it2);
    };

    //!\brief Construct by sequence.
    constexpr aligned_sequence_adaptor_constant_access(inner_type sequence_)
    {
        assign(sequence_);
    };

    //!\brief Construct by initializer_list for inner sequence container.
    constexpr aligned_sequence_adaptor_constant_access(std::initializer_list<value_type> l)
    {
        assign(l);
    };
    //!\}

    /*!\name Iterators
     * \{
    */
    //!\brief Return iterator pointing to first element of underlying sequence.
    auto begin() noexcept
    {
        return iterator{*this, 0};
    }

    //!\brief Return iterator pointing to past-the-end element of gapped sequence.
    auto end() noexcept
    {
        return iterator{*this, size()};
    }

    //!\brief Return const iterator pointing to first element of underlying sequence
    const_iterator cbegin() const noexcept
    {
        return const_iterator{*this, 0};
    }

    //!\brief Return const iterator pointing to past-the-end element of underlying sequence.
    const_iterator cend() const noexcept
    {
        return const_iterator{*this, size()};
    }
     //!\}

    /*!\name Boolean operators
     * \{
     */
    /*!\brief Equality operator for aligned sequences.
    *
    * Two aligned sequences are the same if their literal sequences and gap
    * positions are the same.
    */
    bool operator==(aligned_sequence_t const & rhs) const
    {
        return sequence == rhs.sequence && gap_vector == rhs.gap_vector;
    }

    //!\brief Unequality operator for aligned sequences.
    bool operator!=(aligned_sequence_t const & rhs) const
    {
        return !(*this == rhs);
    }

    //!\brief Swap two aligned sequences and their support structures.
    void swap(aligned_sequence_t & rhs)
    {
        sequence.swap(rhs.sequence);
        std::swap(gap_vector, rhs.gap_vector);
        rs.set_vector(&gap_vector);
        rhs.rs.set_vector(&rhs.gap_vector);
        ss.set_vector(&gap_vector);
        rhs.ss.set_vector(&gap_vector);
    }

    //!\brief Return gapped sequence length.
    size_type size() const
    {
        return gap_vector.size();
    }

    /*!\brief Return the maximal aligned sequence length.
    *
    * The maximal sequence length is limited by either the maximal size of the
    * compressed sequence or the gap vector.
    */
    size_type max_size() const
    {
        return std::min<size_type>(sequence.max_size(), gap_vector.max_size());
    }

    //!\brief An aligned sequence is empty if it contains no alphabet letters.
    bool empty() const
    {
        return sequence.empty();
    }
    //!\}

    /*!\name Assignment
     * \{
    */
    //!\brief Assignment via iterators.
    void assign(iterator it1, iterator it2) // these are aligned_seq iterators!
    {
        sequence.resize(0);
        gap_vector = sdsl::bit_vector(it2 - it1, 0);
        size_type i = 0;
        for (; it1 != it2; ++it1)
        {
            assert(i < this->size());
            if (*it1 == seqan3::gap::GAP)
                gap_vector[i] = 1;
            else
                sequence.push_back(*it1);
            ++i;
        }
        update_support_structures();
    }

    //!\brief Assignment via gapped alphabet sequence.
    void assign(inner_type sequence_)
    {
        sequence.resize(0);
        gap_vector = sdsl::bit_vector(sequence_.size(), 0);
        for (size_type i = 0; i < gap_vector.size(); ++i)
        {
            if (sequence_[i] == gap::GAP)
                gap_vector[i] = 1;
            else
                sequence.push_back(sequence_[i]);
        }
        update_support_structures();
    }

    //!\brief Assignment via initializer_list.
    void assign(std::initializer_list<value_type> l)
    {
        sequence.resize(0);
        gap_vector = sdsl::bit_vector(l.size(), 0);
        for (auto it = std::begin(l); it != std::end(l); ++it){
            if ((*it) == seqan3::gap::GAP)
                gap_vector[it-std::begin(l)] = 1;
            else
                sequence.push_back(*it);
        }
        update_support_structures();
    }

    //!\brief Assign by value repeated 'size' times.
    void assign(size_type size, value_type value)
    {
        if (value == seqan3::gap::GAP){
            gap_vector = sdsl::bit_vector(size, 1);
            sequence.resize(0);
        }
        else {
            gap_vector = sdsl::bit_vector(size, 0);
            sequence.resize(size);
            std::fill(sequence.begin(), sequence.end(), value);
        }
        update_support_structures();
    }
    //!\}

    /*!\name Sequence concept support.
     * \{
    */
    /*!\brief Insert a single value at a given iterator position.
    *
    * Elements right of the newly inserted elemented are shifted right. The
    * returned iterator points to the position of the insertion.
    */
    iterator insert(iterator it, const value_type & value)
    {
        return insert(it, 1, value);
    }

    //!\brief Insert a moved value at a given iterator position.
    iterator insert(iterator it, value_type && value)
    {
        const value_type value_(value);
        return insert(it, value_);
    }

    /*!\brief Insert a value multiple times at an iterator position.
    *
    * The returned iterator points to the position of the left-most inserted
    * element.
    */
    iterator insert(iterator it, size_type size, const value_type & value)
    {
        size_type pos = static_cast<size_type>(it - detail::random_access_iterator<aligned_sequence_t>());
        assert(insert(pos, size, value));
        return it;
    }

    /*!\brief Insert value at a position 'size' times.
    *
    * Return false if given position exceeds the current size by one.
    */
    bool insert(size_type pos, size_type size, const value_type & value)
    {
        if (pos > this->size())
            return false;
        difference_type i, j = static_cast<difference_type>(pos);
        gap_vector.resize(this->size() + size);
        // shift suffix, note that we need i to be a signed integer for the case
        // that the aligned sequence was empty
        for (i = this->size() - size - 1; i >= j; --i)
            gap_vector[i+size] = gap_vector[i];
        // insert gap or alphabet symbols
        if (value == seqan3::gap::GAP)
        {
            for (i = 0; i < static_cast<difference_type>(size); ++i)
                gap_vector[j+i] = 1;
        }
        else{
            for (i = 0; i < static_cast<difference_type>(size); ++i)
                gap_vector[j+i] = 0;
            sequence.insert(sequence.begin() + map_to_underlying_position(j), size, value);
        }
        update_support_structures();
        return true;
    }

    /*\brief Erase element at the iterator's position.
    *
    * Return iterator to past-the-erased element.
    */
    iterator erase(iterator const it)
    {
        assert(erase(static_cast<size_type>(it - iterator{})));
        return it;
    }

    /*!\brief Erase element at a given index.
    *
    * Return false if index exceeds current size minus one.
    */
    bool erase(size_type const pos)
    {
        if (pos >= size())
            return false;
        // erase from compressed sequence string if not gap
        if (!gap_vector[pos])
            sequence.erase(sequence.begin() + map_to_underlying_position(pos));
        // shift suffix left by one and thereby erase i-th bit
        for (size_type i = pos; i < size()-1; ++i)
            gap_vector[i] = gap_vector[i+1];
        gap_vector.resize(size()-1);
        update_support_structures();
        return true;
    }

    //!\brief Erase range between first and second iterator (exclusive).                                              } -> typename type::iterator;
    iterator erase(iterator const it1, iterator const it2)
    {
        size_type i1 = static_cast<size_type>(it1 - iterator{});
        size_type i2 = static_cast<size_type>(it2 - iterator{});
        assert(i1 < this->size() && i2 <= this->size());
        // there is some compressed sequence to be erased, too
        if (rs.rank(i2) - rs.rank(i1) < i2 - i1){
            size_type i1_map = map_to_underlying_position(i1);
            size_type i2_map = map_to_underlying_position(i2-1);
            if (gap_vector[i1]){
                sequence.erase(sequence.begin() + 1 + i1_map, sequence.begin() + 1 + i2_map);
            }else{
                sequence.erase(sequence.begin() + i1_map, sequence.begin() + 1 + i2_map);
            }
        }
        // shift suffix at it2 by the size_type offset it2-it1
        for (size_type i = i2; i < this->size(); ++i){
            gap_vector[i - i2 + i1] = gap_vector[i];
        }
        gap_vector.resize(size() - i2 + i1);
        update_support_structures();
        return it1;
    }

    //!\brief Append an alphabet or gap symbol to the aligned sequence.
    void push_back(value_type symbol)
    {
        assert(max_size() >= size() + 1);
        gap_vector.resize(size() + 1);
        if (symbol == gap::GAP)
            gap_vector[size() - 1] = 1;
        else
            sequence.push_back(symbol);
    }

    //!\brief Delete last element.                                                       } -> void;
    void pop_back()
    {
        assert(this->size() > 0);
        if (!gap_vector[size() - 1])
            sequence.pop_back();
        gap_vector.resize(this->size() - 1);
    }

    //!\brief Clear aligned sequence.
    void clear()
    {
        sequence.clear();
        gap_vector.resize(0);
    }

    //!\brief Return first symbol of aligned sequence.
    value_type front()
    {
        assert(size() > 0u);
        return (*this)[0];
    }

    //!\brief Return last symbol of aligned sequence.
    value_type back()
    {
        assert(size() > 0u);
        return (*this)[size()-1];
    }
    //!\}

    /*!\name Projection between gap and alphabet space.
     * \{
    */
    //!\brief Return gap-free sequence.
    std::vector<value_type> get_underlying_sequence() const
    {
        return sequence;
    }

   /*!\brief Insert gap 'size' times at given index.
   *
   * Return false if given index exceeds sequence length minus one.
   */
   bool insert_gap(size_type pos, size_type size=1)
   {
        if (pos > this->size())
            return false;
        if (pos == this->size() && size == 1)
            push_back(gap::GAP);
        else
            insert(pos, size, gap::GAP);
        return true;
    }

    /*! Remove gap at position pos of length gap_len. Return false when
    *
    * Return false if index exceeds the sequence length or the symbol is not a gap.
    */
    bool remove_gap(size_type pos, size_type len)
    {
        if (pos >= size() || !gap_vector[pos])
            return false;
        gap_vector[pos] = 1;
        return true;
    }

    /*\brief Map a compressed representation index to the aligned sequence index.
    *
    * Note that a i-th 0 corresponds to the i-th position in the compressed sequence.
    * We therefore can directly use the select support structure initialized to
    * count 0s (and not 1s like for the rank) to map to the aligned sequence space.
    * E.g. '--TA-TA--' with input pos = 2 returns select<0>(2) = 5.
    */
    size_type map_to_aligned_position(size_type const pos) const
    {
        assert(pos < gap_vector.size());
        return ss(pos + 1);
    }

    /*!\brief Map from gapped sequence position to index of compressed
    * sequence representation.
    *
    * We use the rank support structure indexed for counting 1s and substract
    * the number gaps in [0; position_gap] from the input position.
    * E.g.           aligned sequence  | - A - - T
    *                 position_gapped  | 0 1 2 3 4
    *      map_to_underlying_position  |-1 0 0 0 1
    * Note that a gap is mapped to the same position than the next preceeding
    * non-gap symbol.
    */
    difference_type map_to_underlying_position(size_type const position_gapped) const
    {
        return static_cast<difference_type>(position_gapped) -
        static_cast<difference_type>(rs.rank(std::min<size_type>(position_gapped+1, this->size())));
    }
    //!\}

    /*!\name Random access sequence concept support.
     * \{
    */
    //!\brief Return reference to aligned sequence for given index.
    constexpr reference operator[](size_type const n) const // const noexcept(noexcept((*host)[pos+n]))
    {
        assert(n < size());
        if (!gap_vector[n]){
            size_type const pos = map_to_underlying_position(n);
            return value_type(sequence[pos]);
        }
        return gap_symbol;
    }

    //!\brief Return reference to aligned sequence for given index.
    value_type at(size_type const idx) const
    {
        return (value_type)(*this)[idx];
    }

    /*!\brief Resize the aligned sequence. Extend with gaps if the new sizes
    * enlarges the current aligned sequence.
    *
    * Note that sdsl does not initialize int_vectors with 0s when applying resize.
    */
    void resize(size_type size)
    {
        size_type size_old = this->size();
        gap_vector.resize(size);
        for (size_type i = size_old; i < this->size(); ++i)
            gap_vector[i] = 1;
        // decrease compressed sequence eventually
        if (size_old > 0)
            sequence.resize(ss.select(size_old-1));
        update_support_structures();
    }
    //{ val.resize(0, typename type::value_type{}) } -> void;
    /*!\brief Resize the aligned sequence. Extend with gaps if the new sizes
    * enlarges the current aligned sequence.
    *
    * Note that sdsl does not initialize int_vectors with 0s when applying resize.
    */
    void resize(size_type size, value_type value)
    {
        if (value == gap::GAP || size < this->size())
            return resize(size);
        // extend with non-gap value size()-size times
        size_type size_old = this->size();
        gap_vector.resize(size);
        for (size_type i = size_old; i < this->size(); ++i)
            gap_vector[i] = 0;
        sequence.insert(sequence.cend(), this->size()-size_old, value);
        update_support_structures();
    }

    //!\}
private:
    //!\privatesection
    //!\brief Where the gap symbol is stored.
    constexpr value_type static const gap_symbol = value_type(gap::GAP);
    /*!\brief Gap presentation of the aligned sequence (1: gap, 0: alphabet letter).
    *
    * The bit vector size corresponds to the true gapped sequence size and is
    * therefore used to answer queries about the gapped sequence size.
    */
    sdsl::bit_vector gap_vector;
    /*!\brief Rank support structure for projection into ungapped sequence space.
    *
    * The rank of position i is number 1s in the prefix 0..i-1. This corresponds
    * to the number of gaps in the aligned sequence. Rank queries are answered in
    * constant time.
    */
    sdsl::rank_support_v5<1,1> rs = sdsl::rank_support_v5<1,1>(&gap_vector);
    /*!\brief Select support structure for projection into gap space.
    *
    * Select i returns the gap vector position of the i-th zero in constant time.
    */
    sdsl::select_support_mcl<0,1> ss = sdsl::select_support_mcl<0,1>(&gap_vector);
    /*!\brief Where the ungapped sequence is stored.
    *
    * Despite the 'compressed' representation of the sequence does not contain
    * gap symbols, its value type is the union the actual alphabet type and the
    * gap symbol. The reason behind is to avoid plenty of castings when assigning
    * or querying for aligned sequence positions.
    */
    std::vector<value_type> sequence{};
};


// global swap { swap(val, val2) } -> void;, TODO: lhs and rhs same inner_type?
template <typename inner_type, char gap_symbol = '_'>
void swap (aligned_sequence_adaptor_constant_access<inner_type> & lhs, aligned_sequence_adaptor_constant_access<inner_type> & rhs)
{
    lhs.swap(rhs);
}


} // namespace seqan3

//static_assert(seqan3::container_concept<seqan3::aligned_sequence_adaptor_constant_access<std::vector<seqan3::gapped_alphabet<seqan3::dna4>>>>);
//static_assert(seqan3::sequence_concept<seqan3::aligned_sequence_adaptor_constant_access<std::vector<seqan3::dna4>>>);
//static_assert(seqan3::random_access_sequence_concept<seqan3::aligned_sequence_adaptor_constant_access<std::vector<seqan3::dna4>>>);
