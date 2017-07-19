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

#include <seqan3/core/concept/core.hpp>
//#include <seqan3/alphabet/alphabet_container.hpp>
#include <seqan3/range/detail/random_access_iterator.hpp>

#include <sdsl/bit_vectors.hpp>
#include "sdsl/rank_support_v5.hpp"
#include "sdsl/select_support_mcl.hpp"

namespace seqan3 {

/*!\brief Implementation of an aligned sequence structure with random access.
 *  \details No iterator operation will modify the container. Arithmetic and boolean
 *  operations are applied to the iterator positions, not the corresponding values
 *  of their containers.
 *  \tparam gapped_alphabet_type The composite alphabet type of the underlying sequence, e.g. dna4 and a gap symbol.
 *  The alphabet_type as part of the gapped_alphabet_t has to satisfy the seqan3::alphabet_concept.
 */

template <template <typename gapped_alphabet_t> struct container_t>
    requires alphabet_concept<gapped_alphabet_t::alphabet_t>
struct aligned_sequence_adaptor_constant_access
{

private:
    //! internal gap representation as bit vector, 0: non-gap, 1: gap
    sdsl::bit_vector gap_vector; // includes support structures for rank/select
    //! support structure to compute select for projection into gap space
    sdsl::select_support_mcl<0> letter_select;
    //! support structure to rank for projection into underlying sequence space
    sdsl::rank_support_v5<1, 1> gap_rank;
    //! TODO: store base sequence without gaps or with gaps?
    using alphabet_t = gapped_alphabet_t::alphabet_t;
    container_t<alphabet_t> sequence;

    // member types required by container_concept
    //!\brief Value type of container elements.
    using value_type = ranges::v3::value_type_t<container_t>;
    //!\brief Use reference type defined by container.
    using reference = typename container_type::reference;
    //!\brief Use const reference type provided by container.
    using const_reference = typename container_type::const_reference;
    //!\brief Use random access iterator on container as iterator type.
    using iterator = random_access_iterator<typename container_type>; //TODO must satisfy forward_iterator_concept and convertible to const_interator
    //!\brief Use const random access iterator on container as const iterator type.
    using const_iterator = random_access_iterator<typename container_type>; //TODO must satisfy forward_iterator_concept
    //!\brief Type for distances between iterators.
    using difference_type = ranges::v3::difference_type_t<container_type>;
    //!\brief Use container's size_type as a position.
    using size_type =  ranges::v3::size_type_t<container_type>;

    //!\brief Use container's size_type as a position.
    // using position_type =  ranges::v3::size_type_t<container_type>;

    // RandomAccessIterator
//    using iterator = detail::ra_iterator<underlying_sequence_type>;
//    using const_iterator = const iterator;
//    using difference_type = uint8_t;

public:

    /*!\name Constructors/Destructors
     * \{
    */
    // \brief Default constructor.
    constexpr aligned_sequence_adaptor_constant_access() = default;

    //!\brief Construct by sequence.
    constexpr aligned_sequence_adaptor_constant_access(container_t & sequence) noexcept : sequence{&sequence} {}

    //!\brief Copy constructor.
    constexpr aligned_sequence_adaptor_constant_access(aligned_sequence_adaptor_constant_access const &) = default;

    //!\brief Copy construction via assignment.
    constexpr aligned_sequence_adaptor_constant_access & operator=(aligned_sequence_adaptor_constant_access const &) = default;

    //!\brief Move constructor.
    constexpr aligned_sequence_adaptor_constant_access (aligned_sequence_adaptor_constant_access &&) = default;

    //!\brief Move assignment.
    constexpr aligned_sequence_adaptor_constant_access & operator=(aligned_sequence_adaptor_constant_access &&) = default;

    //!\brief Use default deconstructor.
    ~aligned_sequence_adaptor_constant_access() = default;
    //!\}
/*
    //! default constructor leaves sequence and gap vector uninitialized
    aligned_sequence_constant_access() :
        gap_vector{sdsl::bit_vector(0, 0)}
        //letter_select{&gap_vector},
        //gap_rank{&gap_vector}
    {
    };

    //optional

    //! init with basic sequence, gap vector remains empty
    aligned_sequence_constant_access(alphabet_container<underlying_sequence_type> sequence_in)
    {
        this->sequence = sequence_in;
        this->gap_vector(bit_vector(0, 0));
        this->letter_select(&this->gap_vector);
        this->gap_rank(&this->gap_vector);
    };
    //! init with basic sequence, gap vector remains empty
    aligned_sequence_constant_access(alphabet_container<underlying_sequence_type> sequence_in,
        rank_support_v5<is_gap, t_pat_len> gap_vector_in) // or bit_vector (see doc of rank_support_v5)
    {
        this->sequence = sequence_in;
        // todo: requirement for bit vector len, or could it be smaller/more basic type, interpretation: no more gaps
        this->gap_vector = gap_vector_in;
        this->letter_select(&this->gap_vector);
        this->gap_rank(&gap_vector);
    };

    //! copy/move constructor
    constexpr aligned_sequence_constant_access(aligned_sequence_constant_access const & rhs) :
        sequence{rhs.sequence},
        gap_vector{rhs.gap_vector}
        //letter_select{&gap_vector},
        //gap_rank{&gap_vector}
    {
    };

    //! assignment operator
    aligned_sequence_constant_access operator=(aligned_sequence_constant_access const & rhs)
    {
        sequence = rhs.sequence;
        gap_vector = rhs.gap_vector;
        //letter_select = &gap_vector;
        //gap_rank = &gap_vector;
        return *this;
    }

    //! desctructor
    ~aligned_sequence_constant_access(void)
    {
        delete sequence;
        delete gap_vector;
        //delete letter_select;
        //delete gap_rank;
    }*/

    /*
    template <typename type>
    concept bool container_concept = requires (type val, type val2)
    {
        // methods and operato
        { val.begin()     } -> typename type::iterator;
        { val.end()       } -> typename type::iterator;
        { val.cbegin()    } -> typename type::const_iterator;
        { val.cend()      } -> typename type::const_iterator;

        { val == val2     } -> bool;
        { val != val2     } -> bool;

        { val.swap(val2)  } -> void;
        { swap(val, val2) } -> void;

        { val.size()      } -> typename type::size_type;
        { val.max_size()  } -> typename type::size_type;
        { val.empty()     } -> bool;
    };
    */

    //! methods and operators required by container concept
    iterator begin()
    {
        return std::begin(sequence);
    }

    iterator end()
    {
        return std::end(sequence);
    }

    const_iterator cbegin() const
    {
        return std::cbegin(sequence);
    }

    const_iterator cend() const
    {
        return std::cend(sequence);
    }

    //! two aligned sequences are the same if their literal sequences and gaps are the same
    bool operator==(aligned_sequence_constant_access const & rhs) const
    {
        return sequence == rhs.sequence && gap_vector == rhs.gap_vector;
    }

    bool operator!=(aligned_sequence_constant_access const & rhs) const
    {
        return sequence != rhs.sequence && gap_vector == rhs.gap_vector;
    }

    /*
    void swap(aligned_sequence_constant_access<alphabet_type> & rhs)
    {
        //sdsl::swap(gap_vector, rhs.gap_vector);
        sequence.swap(rhs.sequence);
        //letter_select.swap(rhs.letter_select);
        //gap_rank.swap(rhs.gap_rank);

    }

    //{ val.max_size()  } -> typename type::size_type;
    //! assume alphabet container as most memory consuming data structure
    size_type max_size()
    {
        return sequence.max_size();
    }

    //{ val.empty()     } -> bool;
    bool empty()
    {
        return sequence.size();
    }

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
template <template <typename template <typename alphabet_t> struct gapped_alphabet_t> struct container_t>
void swap (aligned_sequence_constant_access<gapped_alphabet_t> & lhs, aligned_sequence_constant_access<alphabet_type> & rhs)
{
    lhs.swap(rhs);
}

//static_assert(container_concept<aligned_sequence_constant_access>);
//static_assert(sequence_light_concept<aligned_sequence_constant_access>);

static_assert(static_cast<bool>(container_concept<aligned_sequence_adaptor_constant_access>));
//static_assert(static_cast<bool>(sequence_concept<aligned_sequence_adaptor_constant_access>));
//static_assert(static_cast<bool>(random_access_sequence_concept<aligned_sequence_adaptor_constant_access>));
