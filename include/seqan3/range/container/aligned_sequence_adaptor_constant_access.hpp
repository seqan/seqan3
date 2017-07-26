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

#include <initializer_list>
#include <iostream>
#include <string>
#include <iterator>

#include <seqan3/alphabet/concept.hpp>
#include <seqan3/alphabet/gap/gapped_alphabet.hpp>

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
 *  aligned sequence struct requires the fulfillment of seqan's random access sequence
 *  concept for allowing efficient write and read access.
 *  \tparam gapped_alphabet_type The composite alphabet type of the underlying sequence, e.g. dna4 and a gap symbol.
 *  The alphabet_type as part of the gapped_alphabet_t has to satisfy the seqan3::alphabet_concept.
 *  gap_vector = 1110011000100 and sequence = ATATCGT, gapped sequence is ---AT--ATC-GT
    Note container is define
 */


// todo: require that container is const?

template <typename container_t>
    requires random_access_sequence_concept<container_t> && alphabet_concept<ranges::v3::value_type_t<container_t>>
struct aligned_sequence_adaptor_constant_access
{

private:
    //! internal gap representation as bit vector, 0: non-gap, 1: gap
    sdsl::bit_vector gap_vector{}; // includes support structures for rank/select
    //! support structure to compute select for projection into gap space
    sdsl::select_support_mcl<0> letter_select;
    //! Rank support structure to compute number 1s in prefix 0..i-1.
    sdsl::rank_support_v5<1, 1> rs(&gap_vector);
    //! TODO: store base sequence without gaps or with gaps? dynamic ptr or smart pointer instead? then default constructor ok
    container_t sequence{};

    // specialication of random access iterator
    using gapped_alphabet_t = ranges::v3::value_type_t<container_t>;
    using aligned_sequence_t = aligned_sequence_adaptor_constant_access;
    //detail::iterator_random_access<aligned_sequence> it;
    detail::random_access_iterator<aligned_sequence_t> zero(nullptr, 0);
    detail::random_access_iterator<const aligned_sequence_t> zero_const(nullptr, 0);


public:
    // member types required by container_concept
    //!\brief Value type of container elements.
    using value_type = typename gapped_alphabet_t;
    //!\brief Use reference type defined by container.
    using reference = typename container_t::reference;
    //!\brief Use const reference type provided by container.
    using const_reference = typename container_t::const_reference;
    //!\brief Use random access iterator on container as iterator type.
    // TODO ra iterator over composite data structure, overwrite [], and decrement increment operators
    using iterator = detail::random_access_iterator<aligned_sequence_t>; //TODO must satisfy forward_iterator_concept and convertible to const_interator
    //!\brief Use const random access iterator on container as const iterator type.
    using const_iterator = detail::random_access_iterator<const aligned_sequence_t>; //TODO must satisfy forward_iterator_concept
    //!\brief Type for distances between iterators is taken from alphabet container.
    using difference_type = typename ranges::v3::difference_type_t<container_t>;
    //!\brief Use alphabet container's size_type as a position.
    using size_type = typename ranges::v3::size_type_t<container_t>;

    /*!\name Constructors/Destructors
     * \{
    */
    // \brief Default constructor.
    constexpr aligned_sequence_adaptor_constant_access() = default;

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

    //!\brief Constructors required by sequence concept.
    // TODO: behaviour? replicate input value size times?
    //!\brief Construct by single value repeated size times.
    // TODO feels strange to template with alphabet_t, but construct with gapped_alphabet_t value,
    // add 2nd constructor with alphabet_t as input value?
    constexpr aligned_sequence_adaptor_constant_access(size_type size, value_type value)
    {
        assign(size, value);
    };

    constexpr aligned_sequence_adaptor_constant_access(iterator it1, iterator it2)
    {
        assign(it1, it2);
    };

    //!\brief Construct by sequence.
    //type{std::initializer_list<typename type::value_type>{}}
    constexpr aligned_sequence_adaptor_constant_access(container_t sequence)
    {
        assign(sequence.begin(), sequence.end());
    };

    //{ val = std::initializer_list<typename type::value_type>{}      } -> type &;
    //!\brief Construct by initializer_list for inner sequence container.
    constexpr aligned_sequence_adaptor_constant_access(std::initializer_list<value_type> l)
    {
        assign(l);
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
        // Wrong: would delegate [] to sequence container
        //return iterator(sequence, 0);
        return zero;
    }

    //!\brief Return iterator pointing to past-the-end element of gapped sequence.
    iterator end()
    {
        // TODO: seq.size + rank as 2nd argument, w.r.t. sequence it is outside the alloc mem
        return zero + gap_vector.size();
    }

    //!\brief Return const iterator pointing to first element of underlying sequence
    const_iterator cbegin() const
    {
        return zero_const;
    }

    //!\brief Return const iterator pointing to past-the-end element of underlying sequence.
    const_iterator cend() const
    {
        return zero_const + gap_vector.size();
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

    //!\brief Return gapped sequence length.
    size_type size()
    {
        return sequence.size() + rs.rank(gap_vector.size() - 1;
    }

    //!\brief Maximal aligned sequence size is the one of the ungapped sequence.
    size_type max_size()
    {
        return sequence.max_size(); // plus gap_vector len?
    }

    //!\brief An aligned sequence is empty if the underlying ungapped sequence is empty.
    bool empty()
    {
        return sequence.empty();
    }
    //!\}

    /*!\name Assignment required by sequence concept. replacing its current contents,
     * \{
    */
    //!\brief Assignment via iterators.
    //val.assign(val2.begin(), val2.end())
    void assign(iterator it1, iterator it2) // these are gapped_seq iterators!
    {
        difference_type m = it2 - it1
        sequence.resize(0);
        util::assign(gap_vector, bit_vector(m, 0));
        for (difference_type i = 0; i < m; ++m)
        {
            if (it1[i] == seqan3::gap::gap)
                gap_vector[i] = 1;
            else
                sequence.push_back(it1[i]);
        }
        rs.set_vector(&gap_vector);  // update rank support
    }

    //!\brief Assignment via initializer_list.
    //{ val.assign(std::initializer_list<typename type::value_type>{})
    void assign(std::initializer_list<value_type> l)
    {
        sequence.clear();
        util::assign(gap_vector, bit_vector(l.size(), 0));
        for (auto it = std::begin(l); it != std::end(l); ++it){
            if ((*it) == seqan3::gap::GAP)
                gap_vector[it-std::begin(l)] = 1;
            else
                sequence.push_back(*it);
        }
        rs.set_vector(&gap_vector); // update rank support
    }

    //!\brief Assignment via initializer_list.
    //{ val.assign(typename type::size_type{}, typename type::value_type{}) };
    void assign(size_type size, value_type value)
    {
        if (value == seqan3::gap::GAP)
            util::assign(gap_vector, bit_vector(size, 1));
        else {
            util::assign(gap_vector, bit_vector(size, 0));
            sequence.resize(size);
            std::fill(sequence.begin(), sequence.end(), value);
        }
        rs.set_vector(&gap_vector); // update rank support
    }

    //!\brief Insert single value given by reference at given position. Elements right of insert position are shifted.
    iterator insert(iterator pos, const value_type & value)
    {
        position_type j = pos - zero;
        bit_vector gap_vector_new(gap_vector.size()+1, 0);
        // copy prefix
        for (position_type i = 0; i < j; ++i)
            gap_vector_new[i] = gap_vector[i];

        // insert gap or alphabet symbol
        if (value == seqan3::gap::GAP)
            gap_vector_new[j] = 1;
        else
        {
            position_type j_ungapped = j - rs.rank(j);
            sequence.insert(sequence.begin() + j_ungapped, value);
        }
        // copy suffix
        for (position_type i = j; i < gap_vector.size(); ++i)
            gap_vector_new[i+1] = gap_vector[i];
        // update support structures
        gap_vector.assign(gap_vector_new);
        // TODO: delete gap_vector_new?
        rs.set_vector(&gap_vector); // update rank support
        return pos;
    }

    //{ val.insert(val.begin(), typename type::value_type{})
    iterator insert(iterator pos, value_type && value)
    {
        const value_type value_(std::move(value));
        return insert(pos, value_);
    }

    // { val.insert(val.cbegin(), typename type::size_type{}, typename type::value_type{})} -> typename type::iterator;
    // Insert value multiple times at position indicated by iterator.
    // TODO: more elegant way? Are there any bit shift operators provided by the sdsl?
    iterator insert(const_iterator it1, size_type size, const value_type & value)
    {
        position_type j = pos - zero;
        bit_vector gap_vector_new(gap_vector.size() + size, 0);
        // copy prefix
        position_type i;
        for (i = 0; i < j; ++i)
            gap_vector_new[i] = gap_vector[i];

        // insert gap or alphabet symbols
        if (value == seqan3::gap::GAP)
        {
            for (i = 0; i < size; ++i)
                gap_vector_new[j+i] = 1;
        }
        else
        {
            position_type j_ungapped = j - rs.rank(j);
            sequence.insert(sequence.begin() + j_ungapped + i, size, value);
        }
        // copy suffix
        for (position_type i = j+size; i < gap_vector.size(); ++i)
            gap_vector_new[i+1+size] = gap_vector[i];
        // update support structures
        gap_vector.assign(gap_vector_new);
        rs.set_vector(&gap_vector); // update rank support
        return pos;
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

    //########### random_access_sequence_concept
    //! random access operators
    value_type operator [](size_type idx) const
    {
      // todo: check here for in range or leave it to container?
        if gap_vector[idx]
            return seqan3::gap::GAP;
        difference_type offset = rs.rank(idx);
        return this->sequence[idx - offset];
    }
/*
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
