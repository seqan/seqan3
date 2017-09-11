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
#include <type_traits>
#include <vector>

#include <seqan3/alphabet/concept.hpp>
#include <seqan3/alphabet/gap/gapped.hpp>
#include <seqan3/range/container/concept.hpp>
#include <seqan3/range/detail/random_access_iterator.hpp>

#include <sdsl/sd_vector.hpp>
#include <sdsl/util.hpp>

namespace seqan3 {

/*!\brief Container for an aligned sequence with random access at constant time.
 * \tparam inner_type {The container type with which aligned sequence structure
 * can be constructed, e.g. std::vector<dna4>. The inner type of the
 * container must fulfill the seqan3::alphabet_concept.}
 * \implements partially the seqan3::random_access_sequence_concept
 * \remark {An aligned sequence does not operate on the actual, i.e. ungapped
 * sequence, it holds a pointer to it, and augments it by allowing gap insert or
 * delete operations. Copying the raw initial sequence into the aligned sequence
 * container would be inefficient. The more frequent use case is that a large
 * sequence is loaded once, but remains unmodified in terms of base pair exchanges.
 * Therefore an aligned sequence with random access in constant time only stores
 * a bit vector for the gap information. It is assumed that gaps are rather
 * distributed sparsely. Therefore we use the compressed bit vector sdsl::sd_vector
 * with m*(2+log n/m) space (n = vector length, m = number of set bits).
 * We use rank and a select functions of the sd_vector to map from the gapped
 * sequence space to the ungapped space and vice versa.
 * \par Example
 * '---AT--ATC-GT' is represensed by the bit vector '1110011000100' and a pointer
 * to the ungapped sequence 'ATATCGT'.
 */

// TODO as view, add cp in O(1), bit_vector shared, rename: view_gapped
template <typename inner_type>
//!\cond
    requires alphabet_concept<ranges::v3::value_type_t<inner_type>> &&
    random_access_range_concept<inner_type> && sized_range_concept<inner_type>
//!\endcond
struct aligned_sequence_adaptor_constant_access
{

private:
     //!\privatesection
    using aligned_sequence_t    = aligned_sequence_adaptor_constant_access;
    //using alphabet_t = typename ranges::v3::value_type_t<inner_type>::alphabet_type;
    //!\brief Type of the bit-vector.
    using bit_vector_t          = sdsl::sd_vector<>;
    //!\brief Type of the rank support data structure.
    using rank_support_t        = sdsl::rank_support_sd<1>;     // bit_vector_t::rank_1_type;
    //!\brief Type of the select support data structure.
    using select_support_t      = sdsl::select_support_sd<0>;  //bit_vector_t::select_0_type;

public:
    //!\publicsection
    /*!\name Member types
    * \{
    */
    //!\brief Value type of container elements.
    //!\hideinitializer
    using alphabet_type = typename ranges::v3::value_type_t<inner_type>;
    using value_type = gapped<alphabet_type>;

    //!\brief Use reference type defined by container.
    //!\hideinitializer
    using reference = value_type;

    //!\brief Use const reference type provided by container.
    //!\hideinitializer
    using const_reference = const reference;
    // decltype(std::as_const(data_values)

    //!\brief Use random access iterator on container as iterator type.
    //!\hideinitializer
    using iterator = detail::random_access_iterator<aligned_sequence_t>;

    //!\brief Use const random access iterator on container as const iterator type.
    //!\hideinitializer
    using const_iterator = iterator;

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

    //!\brief Move constructor.
    constexpr aligned_sequence_adaptor_constant_access (aligned_sequence_adaptor_constant_access && rhs) = default;

    //!\brief Move assignment.
    constexpr aligned_sequence_adaptor_constant_access & operator=(aligned_sequence_adaptor_constant_access && rhs) = default;

    //!\brief Use default deconstructor.
    ~aligned_sequence_adaptor_constant_access() = default;
    //!\}

    //!\brief
    /*!\name Constructors of sequence concept.
    * \{
    */
    //!\brief Construct by single value repeated 'size' times.
    constexpr aligned_sequence_adaptor_constant_access(inner_type & sequence) : data{new data_t{&sequence}} {};
//    constexpr random_access_iterator(container_type & host, position_type const pos) noexcept : host{&host}, pos{pos} {}

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
        return data->sequence == rhs.data->sequence && data->gap_vector == rhs.data->gap_vector;
    }

    //!\brief Unequality operator for aligned sequences.
    bool operator!=(aligned_sequence_t const & rhs) const
    {
        return !(*this == rhs);
    }

    //!\brief Swap two aligned sequences and their support structures.
    void swap(aligned_sequence_t & rhs)
    {
        data.swap(rhs.data);
    }

    //!\brief Return gapped sequence length.
    size_type size() const noexcept
    {
        return data->gap_vector.size();
    }

    /*!\brief Return the maximal aligned sequence length.
    *
    * The maximal sequence length is limited by either the maximal size of the
    * compressed sequence or the gap vector. Since the sdsl::sd_vector has no
    * max_size() member function, but can be constructed by a bit_vector, we
    * assume the maximal size is restricted by the one of sdsl::bit_vector.
    */
    size_type max_size() const
    {
        return std::min<size_type>(data->sequence.max_size(), sdsl::bit_vector{}.max_size());
    }

    //!\brief An aligned sequence is empty if it contains no alphabet letters or gaps.
    bool empty() const
    {
        return data->sequence.empty() && data->gap_vector.size() == 0;
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

    /*!\brief Insert a value multiple times at an iterator position.
    *
    * The returned iterator points to the position of the left-most inserted
    * element.
    */
    iterator insert_gap(iterator it, size_type size=1)
    {
        size_type pos = static_cast<size_type>(it - detail::random_access_iterator<aligned_sequence_t>());
        assert(insert_gap(pos, size));
        return it;
    }

    /*!\brief Insert value at a position 'size' times.
    *
    * Return false if given position exceeds the current size by one.
    */
    bool insert_gap(size_type const pos, size_type const size=1)
    {
        if (pos > this->size())
            return false;
        difference_type i, j = static_cast<difference_type>(pos);
        bit_vector_t gap_vector_new = bit_vector_t(sdsl::bit_vector(this->size() + size, 0));
        // copy prefix
        for (i = 0; i < pos; ++i)
            if (data->gap_vector[i])
                gap_vector_new.set(i);

        // insert gap
        for (i = pos; i < pos+size; ++i)
            gap_vector_new.set(i);
        // shift suffix, note that we need i to be a signed integer for the case
        // that the aligned sequence was empty
        for (i = this->size() - size - 1; i >= j; --i)
            if (data->gap_vector[i])
                gap_vector_new.set(i+size);
        data->gap_vector = gap_vector_new;
        // TODO: delete old one?
        data->dirty = true;
        return true;
    }

    /*\brief Erase element at the iterator's position.
    *
    * Return iterator to past-the-erased element.
    */
    iterator erase_gap(iterator const it)
    {
        assert(erase_gap(static_cast<size_type>(it - iterator{})));
        return it;
    }

    /*!\brief Erase element at a given index.
    *
    * Return false if index exceeds current size minus one.
    */
    bool erase_gap(size_type const pos)
    {
        if (pos >= size() || !data->gap_vector[pos])
            return false;
        // init with new size and number of 1 bits
        sd_vector_builder sd_builder{this->size()-1, data->rank_support.rank(this->size()-1)};
        // copy prefix
        for (size_type i = 0; i < pos; ++i)
            if (data->gap_vector[i])
                sd_builder.set(i);
        // shift suffix left by one and thereby erasing i-th bit
        for (size_type i = pos; i < size()-1; ++i)
            if (data->gap_vector[i+1])
                sd_builder.set(i);
        data->gap_vector = bit_vector_t{sd_builder};
        data->dirty = true;
        return true;
    }
    //!\brief Erase all gaps in range pos1 and pos2 (exclusive). Gaps
    // right-hand of 2nd iterator are shifted by the number of gaps deleted.
    bool erase_gap(size_type const pos1, size_type const pos2)
    {
        if (pos1 >= this->size() || pos2 > this->size())
            return false;
        if (data->dirty)
            update_support_structures();
        // number of deleted gaps
        size_type offset = data->select_support.rank(pos2) - data->select_support.rank(pos1);

        sd_vector_builder sd_builder{this->size() - offset, data->rank_support.rank(this->size()) - offset};
        // copy prefix
        for (size_type i = 0; i < pos1; ++i)
            if (data->gap_vector[i])
                sd_builder.set(i);
        // shift suffix at it2 by the size_type offset it2-it1
        for (size_type i = pos1; i < this->size()-offset; ++i)
            if (data->gap_vector[i + offset])
                sd_builder.set(i);
        data->gap_vector = bit_vector_t{sd_builder};
        data->dirty = true;
        return true;
    }

    //!\brief Erase all gaps falling into range given by first and second
    // iterator (exclusive). Gaps right-hand of 2nd iterator are shifted by
    // the number of gaps deleted.
    iterator erase_gap(iterator const it1, iterator const it2)
    {
        size_type pos1 = static_cast<size_type>(it1 - iterator{});
        size_type pos2 = static_cast<size_type>(it2 - iterator{});
        assert(erase_gap(pos1, pos2));
        return it1;
    }

    //!\brief Append gap symbol to the aligned sequence.
    void push_back()
    {
        assert(max_size() >= size() + 1);
        data->gap_vector.resize(size() + 1);
        data->gap_vector[size() - 1] = 1;
        data->dirty = true;
    }

    //!\brief Delete last gap symbol if set, else return false.
    bool pop_back()
    {
        assert(this->size() > 0);
        if (!data->gap_vector[size() - 1])
            return false;
        data->gap_vector.resize(this->size() - 1);
        return true;
    }

    //!\brief Clear gaps in bit vector. Alphabet sequence remains unchanged.
    void clear()
    {
        data->gap_vector = bit_vector_t(sdsl::bit_vector{data->sequence.size(), 0});
        data->dirty = true;
    }

    //!\brief Return first symbol of aligned sequence.
    reference front() const
    {
        assert(size() > 0u);
        return (*this)[0];
    }

    //!\brief Return last symbol of aligned sequence.
    reference back() const
    {
        assert(size() > 0u);
        return (*this)[size()-1];
    }
    //!\}

    /*!\name Projection between gap and alphabet space.
     * \{
    */
    //!\brief Return gap-free sequence.
    inner_type & get_underlying_sequence() const
    {
        return data->sequence;
    }

    void set_underlying_sequence(inner_type & sequence_) const
    {
        data->sequence = sequence_;
        data->gap_vector = bit_vector_t(sdsl::bit_vector{sequence_.size(), 0});
        data->dirty = true;
    }

    /*\brief Map a compressed representation index to the aligned sequence index.
    *
    * Note that a i-th 0 corresponds to the i-th position in the compressed sequence.
    * We therefore can directly use the select support structure initialized to
    * count 0s (and not 1s like for the rank) to map to the aligned sequence space.
    * E.g. '--TA-TA--' with input pos = 2 returns select<0>(2) = 5.
    */
    size_type map_to_aligned_position(size_type const idx)
    {
        //TODO add SEQAN_UNLIKELY
        if (idx >= size())
            throw std::out_of_range{"Trying to access element behind the last in aligned_sequence."};
        if (data->dirty)
            update_support_structures();

        return data->select_support.select(idx + 1);
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
    difference_type map_to_underlying_position(size_type const position_gapped)
    {
        if (data->dirty)
            update_support_structures();

        return static_cast<difference_type>(position_gapped) -
        static_cast<difference_type>(data->rank_support.rank(std::min<size_type>(position_gapped+1, this->size())));
    }
    //!\}

    /*!\name Random access sequence concept support.
     * \{
    */
    //!\brief Return reference to aligned sequence for given index.
    constexpr reference operator[](size_type const idx) // const noexcept(noexcept((*host)[pos+n]))
    {
        assert(idx < size());
        if (!data->gap_vector[idx]){
            size_type const pos = map_to_underlying_position(idx);
            return value_type(data->sequence[pos]);
        }
        return gap::GAP;
    }

    //!\brief Return reference to aligned sequence for given index.
    value_type at(size_type const idx)
    {
        //TODO add SEQAN_UNLIKELY
        if (idx >= size())
            throw std::out_of_range{"Trying to access element behind the last in aligned_sequence."};
        return (value_type)(*this)[idx];
    }

    //!\}
private:
    //!\privatesection
    /*!\brief Gap presentation of the aligned sequence (1: gap, 0: alphabet letter).
    *
    * The bit vector size corresponds to the true gapped sequence size and is
    * therefore used to answer queries about the gapped sequence size.
    */
    struct data_t
    {
        /*!\brief Where the ungapped sequence is stored.
        *
        * Despite the 'compressed' representation of the sequence does not contain
        * gap symbols, its value type is the union the actual alphabet type and the
        * gap symbol. The reason behind is to avoid plenty of castings when assigning
        * or querying for aligned sequence positions. Whether gaps are in the sequence
        * or not will never be tested. Thus gaps in the sequence reference are treated
        * as usual alphabet symbols.
        */
        inner_type & sequence;

        bit_vector_t gap_vector = bit_vector_t();
        /*!\brief Rank support structure for projection into ungapped sequence space.
        *
        * The rank of position i is number 1s in the prefix 0..i-1. This corresponds
        * to the number of gaps in the aligned sequence. Rank queries are answered in
        * constant time.
        */
        //!\hideinitializer
        rank_support_t rank_support{};
        /*!\brief Select support structure for projection into gap space.
        *
        * Select i returns the gap vector position of the i-th zero in constant time.
        */
        //!\hideinitializer
        select_support_t select_support{};

        //!\brief Flag to indicate whether sd_vector support structures needs to be
        // updated before executing rank or select queries.
        bool dirty{true};
    };

    std::shared_ptr<data_t> data;

    //!\brief Re-initialize rank and select support structures gap sd_vector.
    void update_support_structures()
    {
        data->rank_support = sdsl::rank_support_sd<1>(&data->gap_vector);
        data->select_support = sdsl::select_support_sd<0>(&data->gap_vector);
        data->dirty = false;
    }
};

// global swap { swap(val, val2) } -> void;, TODO: lhs and rhs same inner_type?
template <typename inner_type, char gap_symbol = '_'>
void swap (aligned_sequence_adaptor_constant_access<inner_type> & lhs, aligned_sequence_adaptor_constant_access<inner_type> & rhs)
{
    lhs.swap(rhs);
}

} // namespace seqan3
