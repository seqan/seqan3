// ============================================================================
//                 SeqAn - The Library for Sequence Analysis
// ============================================================================
//
// Copyright (c) 2006-2018, Knut Reinert & Freie Universitaet Berlin
// Copyright (c) 2016-2018, Knut Reinert & MPI Molekulare Genetik
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

/*!\file
 * \author Marie Hoffmann <marie.hoffmann AT fu-berlin.de>
 * \brief Contains a gap decorator to annotate sequences with gaps using a set.
 */

#pragma once

#include <iostream>
#include <set>
#include <stdlib.h>

#include <range/v3/all.hpp>

#include <seqan3/alphabet/concept.hpp>
#include <seqan3/alphabet/gap/gapped.hpp>

#define _ 0
#define LOG_LEVEL_AS2 0

// ------------------------------------------------------------------
// gap_decorator_anchor_set
// ------------------------------------------------------------------

namespace seqan3 {

template <typename gap_t>
struct gap_compare {
    bool operator() (const gap_t& lhs, const gap_t& rhs) const {
        return lhs.first < rhs.first;
    }
};

/*!\brief A gap decorator allows the annotation of sequences with gap symbols
 * leaving the underlying sequence unmodified.
 * \implements seqan3::aligned_sequence_concept
 * \ingroup decorator
 *
 * \details
 * The gap_decorator_anchor_set is a modified anchor list approach, i.e. instead
 * of storing anchors and gap lengths relative to the underlying sequence
 * position, the anchor addresses are virtual and gap lengths accumulated from
 * left to right. This reduces the lookup time to log (k) with k being the number
 * of continuous gaps and comes at the price of O(k) worst-case runtime for gap
 * modifications - when inserting or erasing gaps the tailing gaps have to be
 * updated.
 * The gap maintaining structure is a set that is given a gap comparing structure.
 * Sets are implemented as red-black trees and perform random read operations in
 * log (k) time.
 */

// seqan3::aligned_sequence_concept
template <typename inner_type>
struct gap_decorator_anchor_set
{

public:
    /*!\name Member types
    * \{
    */
    using alphabet_type = typename ranges::v3::value_type_t<inner_type>;
    using value_type = gapped<alphabet_type>;
    using reference = value_type;
    using size_type = typename ranges::v3::size_type_t<inner_type>;
    using gap_t = typename std::pair<size_t, size_t>;
    //!\}

    /*!\name Constructors/Destructors
     * \{
     */
    //!\brief Default constructor.
    constexpr gap_decorator_anchor_set()
    {
        data = std::shared_ptr<data_t>(new data_t{});
    };
    //!\brief Copy constructor.
    constexpr gap_decorator_anchor_set(gap_decorator_anchor_set const &) = default;
    //!\brief Copy construction via assignment.
    constexpr gap_decorator_anchor_set & operator=(gap_decorator_anchor_set const &) = default;
    //!\brief Move constructor.
    constexpr gap_decorator_anchor_set (gap_decorator_anchor_set && rhs) = default;
    //!\brief Move assignment.
    constexpr gap_decorator_anchor_set & operator=(gap_decorator_anchor_set && rhs) = default;
    //!\brief Use default deconstructor.
    ~gap_decorator_anchor_set() = default;
    //!\brief Construct by host and explicit position.
    constexpr gap_decorator_anchor_set(inner_type * sequence): data{new data_t{sequence}} {};
    //!\}

    size_type size() const noexcept
    {
        if (data->anchors.size())
            return (*(data->anchors.rbegin())).second + data->sequence->size();
        return data->sequence->size();
    }

    // for benchmark only? delete if needed or allow resizing when tail is exclusively made of gaps
    bool resize(size_type new_size)
    {
        assert(new_size <= this->size());
        for (auto pos = this->size() - 1; pos >= new_size; --pos)
        {
            if ((value_type)(*this)[pos] == gap::GAP)
                erase_gap(pos);
        }
        return true;
    }

    //!\brief Insert a gap of length 'size' at an aligned sequence position 'pos'.
    bool insert_gap(size_type const pos, size_type const size=1)
    {
        typename std::set<gap_t, gap_compare<gap_t>>::iterator it = data->anchors.begin();
        // case 1: extend previous/surrounding gap
        if ((pos < this->size()) && (((value_type)(*this)[pos] == gap::GAP) or (pos > 0 && (value_type)(*this)[pos-1] == gap::GAP)))
        {
            it = data->anchors.lower_bound(gap_t{pos, _});
            if (it == data->anchors.end() || (*it).first > pos)
            {
                it = std::prev(it);
            }
            gap_t gap{(*it).first, (*it).second + size};
            // merge with successor
            auto it_next = it;
            ++it_next;

            if ((*it) < (*std::prev(data->anchors.end())) && (*it_next).first <= (*it).first + size - 1)
            {
                // extend gap for *it, delete *(it+1)
                gap.second += (*it_next).second;
                data->anchors.erase(it_next);
            }
            data->anchors.erase(it);
            data->anchors.insert(gap);

        }
        // case 2: new anchor gap
        else
        {
            gap_t gap{pos, size};
//            data->idx2len[pos] = size;
            // pre: pos not in anchor set, what's the next lower index?
            it = data->anchors.find(gap_t{pos, _}); // return value in case of no lower elem?
            // add accumulated gaps from preceeding gap
            // TODO: same for anchor_set.hpp
            if (it != data->anchors.begin() && it != data->anchors.end()){
                gap.second += (*--it).second;
                //data->idx2len[pos] += data->idx2len[*--it];
            }
            data->anchors.insert(gap);
        }
        // post-processing: reverse update of succeeding gaps
        rupdate(pos, size);
        return true;
    }


    bool erase_gap(size_type const pos)
    {
        return erase_gap(pos, pos+1);
    }

    bool erase_gap(size_type const pos1, size_type const pos2)
    {
        typename std::set<gap_t, gap_compare<gap_t>>::iterator it = data->anchors.lower_bound(gap_t{pos1, _});
        size_type gap_len = get_gap_length(it);

        if (it == data->anchors.end() || (*it).first > pos1)
            it = std::prev(it);
        // case 1: complete gap is deleted
        if (((*it).first == pos1) && (gap_len == pos2-pos1))
            data->anchors.erase(it);
        // case 2: gap to be deleted in tail or larger than 1 (equiv. to shift tail left, i.e. pos remains unchanged)
        else
        {
            gap_t gap{(*it).first, (*it).second - pos2 + pos1};
            // TODO: emplace better?
            data->anchors.erase(it);
            data->anchors.insert(gap);
        }

        // post-processing: forward update of succeeding gaps
        update(pos1, pos2-pos1);
        return true;
    }

    void set_underlying_sequence(inner_type * sequence) const
    {
        data->sequence = sequence;
    }

    constexpr reference operator[](size_type const idx) const
    {
        assert(idx < size());
        // case 1: no gaps
        if (!data->anchors.size()) return value_type((*data->sequence)[idx]);
        // case 2: gaps
        typename std::set<gap_t, gap_compare<gap_t>>::iterator it = data->anchors.lower_bound(gap_t{idx, _});
        if (data->anchors.size() && it == data->anchors.end())
            it = std::prev(it);

        size_type acc = 0, gap_len = 0;
        if ((*it).first <= idx || (it != data->anchors.begin() && (*(std::prev(it))).first <= idx))
        {
            if ((*it).first > idx)    --it;
            acc = (*it).second;
            gap_len = (*it).second;
            if (*it != *(data->anchors.begin()))
                gap_len -= (*(std::prev(it, 1))).second;
            if (idx >= (*it).first && idx < (*it).first + gap_len)
                return gap::GAP;
        }
        return value_type((*data->sequence)[idx - acc]);
    }

private:

    constexpr size_type get_gap_length(typename std::set<gap_t, gap_compare<gap_t>>::iterator it) const noexcept
    {
        if (it == data->anchors.begin()) return (*it).second;
        return (*it).second - (*std::prev(it)).second;
    }

    // reverse update: increase all anchor gaps AFTER position pos by size, i.e. start position AND size
    void rupdate(size_type const pos, size_type const size)
    {
        size_type new_key, new_val;
        for (auto it = std::prev(data->anchors.end(), 1); (*it).first > pos;)
        {
            new_key = (*it).first + size;
            new_val = (*it).second + size;  //data->idx2len[*it] + size;
            data->anchors.insert(gap_t{new_key, new_val});
            data->anchors.erase(*it--);
        }
    }

    //!\brief forward update: decrease all anchor gaps after position pos by size
    void update(size_type const pos, size_type const size)
    {
        assert(pos < size);
        // post: update succeeding gaps  by shifting position key right, start from right to left to avoid collisions
        auto it = data->anchors.lower_bound(gap_t{pos + size + 1, _});

        while (it != data->anchors.end())
        {
            gap_t gap{(*it).first - size, (*it).second - size};
            data->anchors.insert(gap);
            typename std::set<gap_t, gap_compare<gap_t>>::iterator it_next = std::next(it);
            data->anchors.erase(it);
            it = it_next;
        }
    }

    struct data_t
    {
        // pointer to ungapped, underlying sequence
        inner_type * sequence{};
        // store virtual gap positions together with the number of gaps until the given position
        std::set<gap_t, gap_compare<gap_t>> anchors{};
    };
    std::shared_ptr<data_t> data;
};

} // namespace seqan3

#endif
