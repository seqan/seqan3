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
 *  \brief Provides seqn3::gap_decorator_anchor_list.
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

/*!\brief Container for an aligned sequence or to compute an alignment.
 * \par Example
 * '---AT--ATC-GT' is represensed by the bit vector '1110011000100' and a pointer
 * to the ungapped sequence 'ATATCGT'.
 */

template <typename inner_type>
//!\cond
    requires alphabet_concept<ranges::v3::value_type_t<inner_type>> &&
    random_access_range_concept<inner_type> && sized_range_concept<inner_type>
//!\endcond
struct gap_decorator_anchor_list
{

private:
     //!\privatesection
    using gap_decorator_t   = gap_decorator_anchor_list;

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

    using gap_t             = std::pair<size_type, size_type>
    using gap_list_t        = std::vector<gap_t>;

    /* rule of six */
    /*!\name Constructors, destructor and assignment
    * \{
    */
    // \brief Default constructor.
    constexpr gap_decorator_anchor_list()
    {
        inner_type sequence{};
        data = std::shared_ptr<data_t>(new data_t{sequence});
    };

    //!\brief Default copy constructor.
    constexpr gap_decorator_anchor_list(gap_decorator_anchor_list const &) = default;

    //!\brief Default copy construction via assignment.
    constexpr gap_decorator_anchor_list & operator=(gap_decorator_anchor_list const &) = default;

    //!\brief Move constructor.
    constexpr gap_decorator_anchor_list (gap_decorator_anchor_list && rhs) = default;

    //!\brief Move assignment.
    constexpr gap_decorator_anchor_list & operator=(gap_decorator_anchor_list && rhs) = default;

    //!\brief Use default deconstructor.
    ~gap_decorator_anchor_list() = default;
    //!\}

    //!\brief
    /*!\name Constructors of sequence concept.
    * \{
    */
    //!\brief Construct by single value repeated 'size' times
    //shared_ptr<data_t>
    constexpr gap_decorator_anchor_list(inner_type & sequence): data{new data_t{sequence}} {};
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
    * positions are the same. Note there is no operator== in sdsl-lite and the one
    * in the sdsl master branch compares only addresses, but we are interested in
    * the content. This implementation has best case runtime O(1) and worst-case
    * O(r+s+m) where r,s are the initialization costs of sdsl::rank_1_support_sd,
    * and sdsl::select_support_sd, and m the number of 1 bits.
    */
    bool operator==(aligned_sequence_t & rhs) // DONE
    {
        if (data->sequence != rhs.data->sequence || this->size() != rhs.size())
            return false;
        for (std::uint64_t i = 0; i < data->gap_list.size(); ++i)
            if (data->gap_list[i] != rhs.data->gap_list[i])
                return false;
        return true;
    }

    //!\brief Unequality operator for aligned sequences.
    bool operator!=(aligned_sequence_t & rhs)   // DONE
    {
        return !(*this == rhs);
    }

    //!\brief Swap two aligned sequences and their support structures.
    void swap(aligned_sequence_t & rhs)         // DONE
    {
        data.swap(rhs.data);
    }

    //!\brief Return gapped sequence length.
    size_type size() const noexcept             // DONE
    {
        if (!data->gap_list.size())
            return data->sequence.size();
        return data->sequence.size() + data->gap_list.back().second;
    }

    /*!\brief Return the maximal aligned sequence length.
    *
    * The maximal sequence length is limited by either the maximal size of the
    * compressed sequence or the gap vector. Since the sdsl::sd_vector has no
    * max_size() member function, but can be constructed by a bit_vector, we
    * assume the maximal size is restricted by the one of sdsl::bit_vector.
    */
    size_type max_size() const                  // DONE
    {
        return std::min<size_type>(data->sequence.max_size(), data->gap_list.max_size());
    }

    //!\brief An aligned sequence is empty if it contains no alphabet letters or gaps.
    bool empty() const                          // DONE
    {
        return data->sequence.empty() && data->gap_list.size() == 0;
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
    iterator insert_gap(iterator it, size_type size=1)      // DONE
    {
        size_type pos = static_cast<size_type>(it - detail::random_access_iterator<aligned_sequence_t>());
        assert(insert_gap(pos, size));
        return it;
    }

    /*!\brief Insert value at a position 'size' times.
    *
    * Return false if given position exceeds the current size by one.
    */
    bool insert_gap(size_type const pos, size_type const size=1) // TO TEST
    {
        if (pos > this->size())
            return false;
        // case 1: no merging or position computation needed when new gap is first gap
        if (!data->gap_list.size())
            data->gap_list.push_back(gap_t(pos, size));
        else // search true position or expand existing one
        {
            size_type y, x = 0; // current gap range [x .. y[
            bool search_flag = true;
            for (auto it = data->gap_list.begin(); search_flag && it != data->gap_list.end(); ++it)
            {
                x += (*it).first;
                y = x + (*it).second;
                // case 2a: insert before to current gap
                if (pos < x)
                {
                    data->gap_list.insert(it, gap_t(pos, size));
                    search_flag = false;
                }
                // case 2b: pos is gap position or follows directly => expand current gap
                if (pos >= x & pos <= y)
                {
                    *it = gap_t((*it).first, (*it).second + size);
                    search_flag = false;
                }
            }
            // case 2c: new gap starting position is beyond all gaps
            if (search_flag)
                data->gap_list.push_back(gap_t(pos, size));
        }
        return true;
    }

    /* Insert gap relative to sequence, i.e. before the index of the underlying
    * sequence. The maximal allowed insertion position is therefore right after
    * the sequence end. If a gap has already been insert at the given position,
    * false is returned and the state remains unchanged.
    */
    bool insert_gap_rs(size_type const seq_pos, size_type const size=1)
    {
        auto lower = std::lower_bound(data->gap_list.begin(), data->gap_list.end(), seq_pos,
            [seq_pos] (const gap_t& gap) { return seq_pos == gap.first; } );
        // case: push back
        if (lower == data->gap_list.cend())
            data->gap_list.push_back(gap_t(seq_pos, size));
        else if ((*lower).first == seq_pos) // already existing gap
            return false;
        else    // insert before succeeding gap
            data->gap_list.insert(--lower, gap_t(seq_pos, size));
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
    * Worst-case runtime is O(r+s+m). r,s re-initialization of rank and select
    * support, m number of set bits.
    */
    bool erase_gap(size_type const pos)     // DONE
    {
        if ((value_type)(*this)[pos] != gap::GAP)
            return false;
        return erase_gap(pos, pos+1);
    }

    //!\brief Erase all gaps in range pos1 and pos2 (exclusive). Gaps
    /* right-hand of 2nd iterator are shifted by the number of gaps deleted.
    * If the range does not correspond to a contiguous gap, no gaps will be
    * deleted and false will be returned.
    * Remark: this behaviour deviates from the one using a bit_vector, which
    * simply resets all bits in the range without checking connectivity!
    */
    bool erase_gap(size_type const pos1, size_type const pos2)      // UNTESTED
    {
        assert(pos1 <= pos2);
        if (pos1 >= size() || pos2 > size() || !data->gap_list.size() || pos2 < data->gap_list[0].first)
            return false;
        size_type x = 0, y; // current gap range [x; y[
        for (auto it = data->gap_list.begin(); it != data->gap_list.end(); ++it)
        {
            x += (*it).first;
            y = x + (*it).second;
            if (pos1 >= x && pos2 <= y)
            {
                if (pos1 > x) // shorten gap
                    *it = gap_t((*it).first, (*it).second - pos2 + pos1);
                else if (pos1 == x && pos2 == y) // delete gap completely
                    data->gap_list.erase(it);
                else   // remove head of gap
                    *it = gap_t(pos1, (*it).second - pos1 + pos2);
                return true;
            }
        }
        return false;
    }

    //!\brief Erase all gaps falling into range given by first and second
    // iterator (exclusive). Gaps right-hand of 2nd iterator are shifted by
    // the number of gaps deleted.
    iterator erase_gap(iterator const it1, iterator const it2)      // DONE
    {
        size_type pos1 = static_cast<size_type>(it1 - iterator{});
        size_type pos2 = static_cast<size_type>(it2 - iterator{});
        assert(erase_gap(pos1, pos2));
        return it1;
    }

    //!\brief Append gap of length size to the aligned sequence.
    // If last symbol is a gap it will be extended.
    // Remark: pos is global index
    void push_back(size_type pos, size_type size=1)     // UNTESTED
    {
        assert(max_size() >= size() + 1);
        // case 1: there is no gap, sequence is empty or last symbol is not gap
        if (!data->gap_list.size() || !size() || (value_type)(*this)[size()-1] != gap::GAP)
            data->gap_list.push_back(gap_t(pos, size));
        else    // case 2: expand existing one
        {
            auto it = --data->gap_list.end();
            *it = gap_t((*it).first, ++(*it).second);
        }
    }

    // add gap by giving underlying sequence position
    // If a gap has already been inserted at the given position, the existing
    // one will be extended.
    void push_back_rs(size_type pos_rs, size_type size=1)       // UNTESTED
    {
        assert(pos_rs <= data->sequence.size());
        if (!data->gap_list.size())
            data->gap_list.push_back(gap_t(pos_rs, size));
        else if (data->gap_list.back().first == pos_rs)
        {
            auto it = --data->gap_list.end();
            *it = gap_t((*it).first, (*it).second + size);
        }
        else
            data->gap_list.push_back(gap_t());
    }

    //!\brief Delete last gap symbol if set, else return false.
    bool pop_back()
    {
        assert(this->size() > 0);
        if ((value_type)(*this)[this->size()-1] == gap::GAP)
        {   // TODO: one liner by directly writing into the listed pair possible?
            gap_t last = data->gap_list.back();
            data->gap_list[data->gap_list.size()-1] = gap_t(last.first, --last.second);
        }
        else
            return false;
        return true;
    }

    // pop last contiguous gap
    bool pop_back_rs()                  // UNTESTED
    {
        assert(data->gap_list.size() > 0);
        data->gap_list.pop_back();
        return true;
    }

    //!\brief Clear gaps in bit vector. Alphabet sequence remains unchanged.
    void clear()                        // UNTESTED
    {
        data->gap_list.clear();
    }

    //!\brief Return first symbol of aligned sequence.
    reference front()                                           // DONE
    {
        assert(size() > 0u);
        return (*this)[0];
    }

    //!\brief Return last symbol of aligned sequence.
    reference back()                                            // DONE
    {
        assert(size() > 0u);
        return (*this)[size()-1];
    }
    //!\}

    /*!\name Sequence getter and setter.
     * \{
    */
    //!\brief Return pointer to gap-free sequence.
    inner_type & get_underlying_sequence() const                // DONE
    {
        return data->sequence;
    }

    //!\brief Set pointer to ungapped sequence and reset gap vector.
    void set_underlying_sequence(inner_type & sequence) const   // DONE
    {
        data->sequence = sequence;
        data->gap_list.clear();
    }
    //!\}

    /*\brief Map a compressed representation index to the aligned sequence index.
    *
    * Note that a i-th 0 corresponds to the i-th position in the compressed sequence.
    * E.g. '--TA-TA--' with input pos = 2 returns 5.
    */
    size_type map_to_aligned_position(size_type const idx)      // UNTESTED
    {
        //TODO add SEQAN_UNLIKELY
        if (idx >= data->sequence.size())
            throw std::out_of_range{"Trying to access element behind the last in aligned_sequence."};
        // case 1: no gap before position idx
        if (!data->gap_list.size() || data->gap_list[0].first < idx)
            return idx;
        // case 2: sum up gaps until idx-th underlying sequence position
        size_type offset = 0; // gap sum up to position idx in anchor gap vector
        auto upper = std::upper_bound(data->gap_vector.begin(), data->gap_vector.end(),
            idx, [idx] (gap_t gap) {return idx < gap.first});
        size_type s = std::accumulate(data->gap_vector.begin(), upper);
        return s;
    }

    /*!\name Random access sequence concept support.
     * \{
    */
    //!\brief Return reference to aligned sequence for given index.
    constexpr reference operator[](size_type const idx) // const noexcept(noexcept((*host)[pos+n]))
    {
        assert(idx < size());
        if (!data->gap_list.size() || position_gapped < data->gap_list[0].first)
            return position_gapped;
        // case 2: compute position in gap or between two gaps
        auto it = data->gap_list.begin();
        difference_type acc = (*it).first;
        ++it;
        // sum up gap offsets and gap lengths
        while (it != data->gap_list.end() && idx > acc + (*it).first - (*(it-1)).first + (*it).second)
        {
            acc += (*it).first - (*(it-1)).first + (*it).second;
            ++it;
        }
        if (idx >= acc + (*it).first - (*(it-1)).first)
            return gap::GAP;
        else
            return data->sequence[idx - acc];
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
    //!\brief Structure for storing a sequence and gap information and helper functions.

    struct data_t
    {
        /*!\brief Where the ungapped sequence is stored.
        *
        * The ungapped sequence is the original sequence of arbitrary alphabet type.
        * If the alphabet type allows gap symbols, these are treated as normal symbols.
        * Only gaps inserted via this interface are stored in a bit vector.
        */
        inner_type & sequence;

        /*!\brief Where the gaps are stored.
        *
        */
        gap_list_t gap_list = gap_list_t();

    };
    std::shared_ptr<data_t> data;
};

//!\brief Global swap function.
template <typename inner_type, char gap_symbol = '_'>
void swap (gap_decorator_anchor_list<inner_type> & lhs, gap_decorator_anchor_list<inner_type> & rhs)
{
    lhs.swap(rhs);
}

} // namespace seqan3
