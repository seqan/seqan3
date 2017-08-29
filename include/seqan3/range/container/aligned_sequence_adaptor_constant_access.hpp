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

    #include <algorithm>
    #include <initializer_list>
    #include <iostream>
    #include <string>
    #include <iterator>

#include <seqan3/alphabet/concept.hpp>
#include <seqan3/alphabet/gap/gapped.hpp>

//#include <seqan3/core/concept/core.hpp>
//#include <seqan3/alphabet/alphabet_container.hpp>
#include <seqan3/range/container/concept.hpp>
#include <seqan3/range/detail/random_access_iterator.hpp>

#include <sdsl/bit_vectors.hpp>
//#include <sdsl/bit_vectors.hpp>

#include "sdsl/rank_support_v5.hpp"
#include "sdsl/select_support_mcl.hpp"
#include "sdsl/util.hpp"

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

template <typename container_t>
    //requires random_access_sequence_concept<container_t> && alphabet_concept<ranges::v3::value_type_t<container_t>>
    requires alphabet_concept<ranges::v3::value_type_t<container_t>> // TODO: think ra on container needed?
struct aligned_sequence_adaptor_constant_access
{

private:
    using aligned_sequence_t = aligned_sequence_adaptor_constant_access;
    //using alphabet_t = typename ranges::v3::value_type_t<container_t>::alphabet_type;

// re-initialize rank support structure of bit_vector
void update_support_structures()
{
    rs = sdsl::rank_support_v5<1, 1>(&gap_vector);
    ss = sdsl::select_support_mcl<0, 1>(&gap_vector);
    //rs.set_vector(&gap_vector);
}

public:
    // member types required by container_concept
    //!\brief Value type of container elements.
    using value_type = typename ranges::v3::value_type_t<container_t>;
    //!\brief Use reference type defined by container.
    using reference = value_type;
    //!\brief Use const reference type provided by container.
    using const_reference = const reference;
    //!\brief Use random access iterator on container as iterator type.
    using iterator = detail::random_access_iterator<aligned_sequence_t>;
    //!\brief Use const random access iterator on container as const iterator type.
    using const_iterator = detail::random_access_iterator<aligned_sequence_t const>;
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
    constexpr aligned_sequence_adaptor_constant_access (aligned_sequence_adaptor_constant_access && rhs) :
    gap_vector(std::move(rhs.gap_vector)), sequence(std::move(rhs.sequence))
    {
        rs.set_vector(&gap_vector);
    } // why default implicitly deleted?

    //!\brief Move assignment.
    constexpr aligned_sequence_adaptor_constant_access & operator=(aligned_sequence_adaptor_constant_access && rhs)
    {
        gap_vector = std::move(rhs.gap_vector);
        sequence = std::move(rhs.sequence);
        rs.set_vector(&gap_vector);
        return *this;
    } //= default;

    //!\brief Use default deconstructor.
    ~aligned_sequence_adaptor_constant_access() = default;

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
    constexpr aligned_sequence_adaptor_constant_access(container_t sequence_)
    {
        assign(sequence_);
    };

    //{ val = std::initializer_list<typename type::value_type>{}      } -> type &;
    //!\brief Construct by initializer_list for inner sequence container.
    constexpr aligned_sequence_adaptor_constant_access(std::initializer_list<value_type> l)
    {
        assign(l);
    };

    //!\}

    /*!\name Container concept requirements.
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

    //! two aligned sequences are the same if their literal sequences and gaps are the same
    bool operator==(aligned_sequence_t const & rhs) const
    {
        return sequence == rhs.sequence && gap_vector == rhs.gap_vector;
    }

    bool operator!=(aligned_sequence_t const & rhs) const
    {
        return !(*this == rhs);
    }

    //!\brief Swap all relevant sequence structures.
    // Relevant are gap-free sequence, bit vector, and rank support
    void swap(aligned_sequence_t & rhs)
    {
        sequence.swap(rhs.sequence);
        std::swap(gap_vector, rhs.gap_vector);
//        std::swap(rs, rhs.rs);
        // hand over new addresses of bit_vector to rank support structure explicitly
        // more: https://github.com/xxsds/sdsl-lite/issues/15
        rs.set_vector(&gap_vector);
        rhs.rs.set_vector(&rhs.gap_vector);
    }

    //!\brief Return gapped sequence length.
    size_type size() const
    {
        return gap_vector.size();
    }

    //!\brief Maximal aligned sequence size is the one of the ungapped sequence.
    size_type max_size() const
    {
        // max_size equal to one of ungapped sequence
        return std::min<size_type>(sequence.max_size(), gap_vector.max_size()); // plus gap_vector len?
    }

    //!\brief An aligned sequence is empty if the underlying ungapped sequence is empty.
    bool empty() const
    {
        return sequence.empty();
    }

    //!\}

    /*!\name Assignment required by sequence concept. replacing its current contents,
     * \{
    */
    //!\brief Assignment via iterators.
    //val.assign(val2.begin(), val2.end())

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

    //!\brief Assignment via gapped alphabet sequence. Note this is not required
    // by the container concept.
    void assign(container_t sequence_)
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
    //{ val.assign(std::initializer_list<typename type::value_type>{})
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

    //!\brief Assignment via initializer_list.
    //{ val.assign(typename type::size_type{}, typename type::value_type{}) };
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


    //!\brief Insert single value given by reference at given position.
    // Elements right of insert position are shifted.
    iterator insert(iterator it, const value_type & value)
    {
        return insert(it, 1, value);
    }

    //{ val.insert(val.begin(), typename type::value_type{})
    iterator insert(iterator it, value_type && value)
    {
        const value_type value_(value);
        return insert(it, value_);
    }


    // { val.insert(val.cbegin(), typename type::size_type{}, typename type::value_type{})} -> typename type::iterator;
    // Insert value multiple times at position indicated by iterator.
    // TODO: more elegant way? Are there any bit shift operators provided by the sdsl?
    // return iterator to left most inserted element
    iterator insert(iterator it, size_type size, const value_type & value)
    {
        size_type pos = static_cast<size_type>(it - detail::random_access_iterator<aligned_sequence_t>());
        assert(insert(pos, size, value));
        return it;
    }

    bool insert(size_type pos, size_type size, const value_type & value)
    {
        if (pos > this->size())
            return false;
        difference_type i, j = static_cast<difference_type>(pos);
        difference_type size_old = this->size();
        gap_vector.resize(this->size() + size);
        // shift suffix, note we need i to be a signed integer
        for (i = size_old - 1; i >= j; --i)
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
        // update support structures
        update_support_structures();
        return true;
    }


    // Erase element at the iterator's position and return an iterator pointing
    // to the new location of the element that followed the last element erased.
    //{ val.erase(val.cbegin())                   } -> typename type::iterator;
    iterator erase(iterator const it)
    {
        size_type pos = static_cast<size_type>(it - iterator{});
        assert(erase(pos));
        return it;
    }

    bool erase(size_type const pos)
    {
        if (pos >= size())
            return false;
        // erase from compressed sequence string
        if (!gap_vector[pos])
            sequence.erase(sequence.begin() + map_to_underlying_position(pos));
        // shift suffix left
        for (size_type i = pos; i < size()-1; ++i)
            gap_vector[i] = gap_vector[i+1];
        // update rank support only there are 1s after erased position
        if (rs.rank(pos+1) != rs.rank(size()))
            rs.set_vector(&gap_vector);
        gap_vector.resize(size()-1);
        return true;
    }


//    { val.erase(val.cbegin(), val.cend())                                              } -> typename type::iterator;
    iterator erase(iterator const pos1, iterator const pos2)
    {
        size_type i1 = static_cast<size_type>(pos1 - iterator{});
        size_type i2 = static_cast<size_type>(pos2 - iterator{});
        assert(i1 < this->size() && i2 <= this->size());
        if (rs.rank(i2) - rs.rank(i1) < i2 - i1){  // there is some compressed sequence to be erased, too
            size_type i1_map = map_to_underlying_position(i1);
            size_type i2_map = map_to_underlying_position(i2-1);
            if (gap_vector[i1]){
                sequence.erase(sequence.begin() + 1 + i1_map, sequence.begin() + 1 + i2_map);
            }else{
                sequence.erase(sequence.begin() + i1_map, sequence.begin() + 1 + i2_map);
            }
        }
        // copy gap bits if there are bits set beyond pos1
        for (size_type i = i2; i < this->size(); ++i){
            gap_vector[i - i2 + i1] = gap_vector[i];
        }
        gap_vector.resize(size() - i2 + i1);
        update_support_structures(); //rs.set_vector(&gap_vector);
        return pos1;
    }

//    { val.push_back(val.front())                                                       } -> void;
    void push_back(value_type symbol)
    {
        gap_vector.resize(size() + 1);
        if (symbol == gap::GAP)
            gap_vector[size() - 1] = 1;
        else
            sequence.push_back(symbol);
    }

    // same as above?{ val.push_back(typename type::value_type{})                                       } -> void;
    //{ val.pop_back()                                                                   } -> void;
    void pop_back()
    {
        assert(this->size() > 0);
        if (!gap_vector[size() - 1])
            sequence.pop_back();
        gap_vector.resize(this->size() - 1);
    }

    //{ val.clear()
    void clear()
    {
        sequence.clear();
        gap_vector.resize(0);
    }

    //{ val.front() } -> typename type::value_type &;
    value_type front()
    {
        return (*this)[0];
    }

    //    { val.back()  } -> typename type::value_type &;
    value_type back()
    {
        return (*this)[size()-1];
    }
    //!\}

    /*!\name Aligned sequence functions for moving between gapped and alphabet space.
     * \{
    */
    //! return gap-free sequence
    std::vector<value_type> get_underlying_sequence() const
    {
        return sequence;
    }

   //! insert a gap at position pos of length gap_len
   // insertion at position of the current size is a push_back
   // overload insert?, because insert(iterator, ...) transforms back into pos
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

    //! Remove gap at position pos of length gap_len. Return false when
    // to be removed gap positions exceed sequence length or are not gap,
    // keep sequence then unchanged, else remove and return true.
    //
    bool remove_gap(size_type pos, size_type len)
    {
        if (pos >= size() || !gap_vector[pos])
            return false;
        gap_vector[pos] = 1;
        return true;
    }

    // rank computation to compute projection from gap to sequence space
    //! Given the underlying sequence position project into gap
    // space by adding the gap rank.
    // '--TA-TA--' with position_base=5 returns 5-rnk(5) = 2
    size_type map_to_aligned_position(size_type const pos) const
    {
        //! return index after #pos aligned letters, boundary check by sdsl
        assert(pos < gap_vector.size());
        return ss(pos + 1);
    }


    //! Given the aligned sequence position (possibly including gaps) project
    // into gap-free alphabet space by removing the gap rank.
    // substract number of gaps in [0; position_gap], if given position is a gap,
    // e.g.           aligned sequence  | - A - - T
    //                 position_gapped  | 0 1 2 3 4
    //      map_to_underlying_position  |-1 0 0 0 1
    // Mapping gaps to the compressed sequence representation is ambiguous and therefore asserted
    difference_type map_to_underlying_position(size_type const position_gapped) const
    {
        return static_cast<difference_type>(position_gapped) - static_cast<difference_type>(rs.rank(std::min<size_type>(position_gapped+1, this->size())));
    }
    //!\}

    //########### random_access_sequence_concept
    //! random access operators
    // TODO: check here for in range or leave it to container?
    // TODO: is this still a constexpr, there is a function call which may not be constexpr
    constexpr reference operator[](size_type const n) const // const noexcept(noexcept((*host)[pos+n]))
    {
        assert(n < size());
        if (!gap_vector[n]){
            size_type const pos = map_to_underlying_position(n);
            return value_type(sequence[pos]);
        }
        return gap_symbol;
    }

    //! random access operators
    value_type at(size_type const idx) const
    {
        return (value_type)(*this)[idx];
    }

    // modify container
    //{ val.resize(0)                              } -> void;
    void resize(size_type size)
    {
        // convert to ungapped sequence index
        // reset all bits beyond position size
        //for (difference_type)
        //sequence.resize(size);
    }
    //{ val.resize(0, typename type::value_type{}) } -> void;

private:
    // TODO: if gap symbol for 2 instance can be different, then provide get_gap fct
    constexpr value_type static const gap_symbol = value_type(gap::GAP); // seqan3::gap::GAP
    //size_type gapped_sequence_size{0};
    //!\brief gap_vector stores for each virtual gapped sequence position either a 0 (non-gap) or a 1 (gap).
    // its size corresponds therefore to the true gapped sequence size. Whereas sequence is the packed letter series.
    sdsl::bit_vector gap_vector;// = sdsl::bit_vector(bit_vector_size, 0); // includes support structures for rank/select

    //! Rank support structure to compute number 1s in prefix 0..i-1.
    sdsl::rank_support_v5<1,1> rs = sdsl::rank_support_v5<1,1>(&gap_vector);
    //! support structure to compute select for projection into gap space
    //template <uint8_t t_bit_pattern, uint8_t t_pattern_len>
    sdsl::select_support_mcl<0,1> ss = sdsl::select_support_mcl<0,1>(&gap_vector);
    //sdsl::bit_vector::select_0_type ss(&gap_vector);

    //! TODO: store base sequence without gaps or with gaps? dynamic ptr or smart pointer instead? then default constructor ok
    //!\brief: ungapped sequence.
    // Despite the 'compressed' representation of the sequence does not contain gap symbols,
    // it is vector of type value_type which is the union of an alphabet type and the gap symbol type.
    // The reason why it is not a pure alphabet type is to avoid plenty of useless castings
    // between gapped and raw alphabet types and function overloads.
    std::vector<value_type> sequence{};
};


// global swap { swap(val, val2) } -> void;, TODO: lhs and rhs same container_t?
template <typename container_t, char gap_symbol = '_'>
void swap (aligned_sequence_adaptor_constant_access<container_t> & lhs, aligned_sequence_adaptor_constant_access<container_t> & rhs)
{
    lhs.swap(rhs);
}


} // namespace seqan3

//static_assert(seqan3::container_concept<seqan3::aligned_sequence_adaptor_constant_access<std::vector<seqan3::gapped_alphabet<seqan3::dna4>>>>);
//static_assert(seqan3::sequence_concept<seqan3::aligned_sequence_adaptor_constant_access<std::vector<seqan3::dna4>>>);
//static_assert(seqan3::random_access_sequence_concept<seqan3::aligned_sequence_adaptor_constant_access<std::vector<seqan3::dna4>>>);
