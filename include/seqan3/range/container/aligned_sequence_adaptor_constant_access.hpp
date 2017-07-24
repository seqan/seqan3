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

/*! \file container/aligned_sequences_constant_access.hpp
 *  \brief Aligned sequences with random access operator.
 *  \ingroup container
 *  \author Marie Hoffmann <marie.hoffmann AT fu-berlin.de>
 */

#pragma once

#include <iostream>
#include <string>
#include <iterator>

#include <seqan3/alphabet/concept.hpp>

//#include <seqan3/core/concept/core.hpp>
//#include <seqan3/alphabet/alphabet_container.hpp>
#include <seqan3/range/container/concept.hpp>
#include <seqan3/range/detail/random_access_iterator.hpp>

#include <sdsl/bit_vectors.hpp>
#include "sdsl/rank_support_v5.hpp"
#include "sdsl/select_support_mcl.hpp"

namespace seqan3 {

/*!\brief Implementation of an aligned sequence structure with random access.
 *  \details No iterator operation will modify the container. Arithmetic and boolean
 *  operations are applied to the iterator positions, not the corresponding values
 *  of their containers. The aligned_sequence with constant access operator has to fulfill the
 *  seqan3::aligned_sequence_adaptor_constant_access concept which includes the
 *  seqan's container and sequence concepts. The ungapped sequence stored by the
 *  aligned sequence struct requires only seqan's sequence concept for assignment and clearing operations.
 *  \tparam gapped_alphabet_type The composite alphabet type of the underlying sequence, e.g. dna4 and a gap symbol.
 *  The alphabet_type as part of the gapped_alphabet_t has to satisfy the seqan3::alphabet_concept.
 */


// todo: require that container is const?
template <typename container_t, char gap_symbol='_'>
    requires sequence_concept<container_t> && alphabet_concept<ranges::v3::value_type_t<container_t>>
struct aligned_sequence_adaptor_constant_access
{

private:
    //! internal gap representation as bit vector, 0: non-gap, 1: gap
    sdsl::bit_vector gap_vector{}; // includes support structures for rank/select
    //! support structure to compute select for projection into gap space
    sdsl::select_support_mcl<0> letter_select;
    //! support structure to rank for projection into underlying sequence space
    sdsl::rank_support_v5<1, 1> gap_rank;
    //! TODO: store base sequence without gaps or with gaps? dynamic ptr or smart pointer instead? then default constructor ok
    container_t sequence{};

    // Shortcuts
    using aligned_sequence = aligned_sequence_adaptor_constant_access;

public:
    // member types required by container_concept
    //!\brief Value type of container elements.
    using value_type = typename ranges::v3::value_type_t<container_t>;
    //!\brief Use reference type defined by container.
    using reference = typename container_t::reference;
    //!\brief Use const reference type provided by container.
    using const_reference = typename container_t::const_reference;
    //!\brief Use random access iterator on container as iterator type.
    using iterator = detail::random_access_iterator<container_t>; //TODO must satisfy forward_iterator_concept and convertible to const_interator
    //!\brief Use const random access iterator on container as const iterator type.
    using const_iterator = detail::random_access_iterator<const container_t>; //TODO must satisfy forward_iterator_concept
    //!\brief Type for distances between iterators.
    using difference_type = typename ranges::v3::difference_type_t<container_t>;
    //!\brief Use container's size_type as a position.
    using size_type = typename ranges::v3::size_type_t<container_t>;

    /*!\name Constructors/Destructors
     * \{
    */
    // \brief Default constructor.
    constexpr aligned_sequence_adaptor_constant_access() = default;

    //!\brief Construct by sequence.
    // TODO: init helpers
    constexpr aligned_sequence_adaptor_constant_access(container_t sequence) noexcept : sequence{sequence} {};

    //!\brief Copy constructor.
    constexpr aligned_sequence_adaptor_constant_access(aligned_sequence_adaptor_constant_access const &) = default;

    //!\brief Copy construction via assignment.
    constexpr aligned_sequence_adaptor_constant_access & operator=(aligned_sequence_adaptor_constant_access const &) = default;

    //!\brief Move constructor.
    constexpr aligned_sequence_adaptor_constant_access (aligned_sequence_adaptor_constant_access &&) = default;

    //!\brief Move assignment.
    constexpr aligned_sequence_adaptor_constant_access & operator=(aligned_sequence_adaptor_constant_access &&) = default;

    //!\brief Use default deconstructor.
    ~aligned_sequence_adaptor_constant_access()
    {
            sequence.clear();
            sequence.shrink_to_fit();
            // TODO: clear other data structures
    }; //= default;

    //!\brief Constructor required by sequence concept.
// { type{typename type::size_type{}, typename type::value_type{}} };
    // TODO: behaviour? replicate input value size times?
    constexpr aligned_sequence_adaptor_constant_access(size_type size, value_type value)
    {
        container_t sequence;
    };


    //!\}
/*
    //! default constructor leaves sequence and gap vector uninitialized
    aligned_sequence_adaptor_constant_access() :
        gap_vector{sdsl::bit_vector(0, 0)}
        //letter_select{&gap_vector},
        //gap_rank{&gap_vector}
    {
    };

    //optional

    //! init with basic sequence, gap vector remains empty
    aligned_sequence(alphabet_container<underlying_sequence_type> sequence_in)
    {
        this->sequence = sequence_in;
        this->gap_vector(bit_vector(0, 0));
        this->letter_select(&this->gap_vector);
        this->gap_rank(&this->gap_vector);
    };
    //! init with basic sequence, gap vector remains empty
    aligned_sequence(alphabet_container<underlying_sequence_type> sequence_in,
        rank_support_v5<is_gap, t_pat_len> gap_vector_in) // or bit_vector (see doc of rank_support_v5)
    {
        this->sequence = sequence_in;
        // todo: requirement for bit vector len, or could it be smaller/more basic type, interpretation: no more gaps
        this->gap_vector = gap_vector_in;
        this->letter_select(&this->gap_vector);
        this->gap_rank(&gap_vector);
    };
*/

    /*!\name Container concept requirements.
     * \{
    */
    //!\brief Return iterator pointing to first element of underlying sequence.
    iterator begin()
    {
        return iterator(sequence, 0); //std::begin(sequence);
    }

    //!\brief Return iterator pointing to past-the-end element of underlying sequence.
    iterator end()
    {
        return iterator(sequence, sequence.size());
    }

    //!\brief Return const iterator pointing to first element of underlying sequence
    const_iterator cbegin() const
    {
        return const_iterator(sequence, 0);;
    }

    //!\brief Return const iterator pointing to past-the-end element of underlying sequence.
    const_iterator cend() const
    {
        return const_iterator(sequence, sequence.size());
    }

    //! two aligned sequences are the same if their literal sequences and gaps are the same
    bool operator==(aligned_sequence const & rhs) const
    {
        return sequence == rhs.sequence && gap_vector == rhs.gap_vector;
    }

    bool operator!=(aligned_sequence const & rhs) const
    {
        return sequence != rhs.sequence && gap_vector == rhs.gap_vector;
    }

    void swap(aligned_sequence & rhs)
    {
        //sdsl::swap(gap_vector, rhs.gap_vector);
        sequence.swap(rhs.sequence);
        // TODO: swap gap data structures
        //letter_select.swap(rhs.letter_select);
        //gap_rank.swap(rhs.gap_rank);
    }

    //!\brief Return ungapped sequence length.
    size_type size()
    {
        return sequence.size();
    }

    //!\brief Maximal aligned sequence size is the one of the ungapped sequence.
    size_type max_size()
    {
        return sequence.max_size();
    }

    //!\brief An aligned sequence is empty if the underlying ungapped sequence is empty.
    bool empty()
    {
        return sequence.empty();
    }
    //!\}

/*
    //! return gap-free sequence
    vector<alphabet_type> get_underlying_sequence() const
    {
        return sequence;
    }

    //! insert a gap at position pos of length gap_len
    bool insert_gap(size_type pos, size_type len)
    {

    }

    //! Remove gap at position pos of length gap_len. Return false when
    // to be removed gap positions exceed sequence length or are not gap,
    // keep sequence then unchanged, else remove and return true.
    //
    bool remove_gap(size_type pos, size_type len)
    {

    }

    // rank computation to compute projection from gap to sequence space
    //! Given the underlying sequence position project into gap
    // space by adding the gap rank.
    //
    size_type map_to_aligned_position(size_type position_base)
    {
        //! return index after #position_base aligned letters
        return letter_select(position_base);
    }

    //! Given the aligned sequence position (possibly including gaps) project
    // into gap-free alphabet space by removing the gap rank.
    //
    size_type map_to_underlying_position(size_type position_gap)
    {
        // substract number of gaps in [0; position_gap]
        return position_gap - gap_rank();
    }
*/

    // access container
    //{ val[0]    } -> typename type::value_type &;
    //{ val.at(0) } -> typename type::value_type &;

/*

    //########### random_access_sequence_concept
    //! random access operators
    value_type operator [](size_type idx) const
    {
      // todo: check here for in range or leave it to container?
        return this->sequence[idx];
    }

    //! random access operators
    value_type at(size_type idx) const
    {
      // todo: check here for in range or leave it to container?
        return this->sequence[idx];
    }

    // modify container
    void resize(size_type)                              } -> void;
    { val.resize(0, typename type::value_type{}) } -> void;
*/
};

// global swap { swap(val, val2) } -> void;, TODO: lhs and rhs same container_t?
template <typename container_t, char gap_sybmol = '_'>
void swap (aligned_sequence_adaptor_constant_access<container_t> & lhs, aligned_sequence_adaptor_constant_access<container_t> & rhs)
{
    lhs.swap(rhs);
    // TODO: swap other relevant data structures
}

} // namespace seqan3

static_assert(seqan3::container_concept<seqan3::aligned_sequence_adaptor_constant_access<std::vector<seqan3::dna4>>>);
//static_assert(seqan3::sequence_concept<seqan3::aligned_sequence_adaptor_constant_access<std::vector<seqan3::dna4>>>);
//static_assert(seqan3::random_access_sequence_concept<seqan3::aligned_sequence_adaptor_constant_access<std::vector<seqan3::dna4>>>);
