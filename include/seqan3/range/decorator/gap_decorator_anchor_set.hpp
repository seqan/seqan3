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
  * \brief Contains gap decorator to annotate sequences with gaps using a set.
  */

#pragma once

#include <set>
#include <tuple>
#include <type_traits>

#include <range/v3/all.hpp>

#include <seqan3/alignment/aligned_sequence/aligned_sequence_concept.hpp>
#include <seqan3/alphabet/concept.hpp>
#include <seqan3/alphabet/gap/gapped.hpp>
#include <seqan3/range/container/concept.hpp>
#include <seqan3/range/detail/random_access_iterator.hpp>

namespace seqan3
{

/*!\brief A gap decorator allows the annotation of sequences with gap symbols
 * leaving the underlying sequence unmodified.
 * \tparam inner_type The type of alphabet sequences that will be stored.
 * \implements seqan3::aligned_sequence_concept
 * \ingroup decorator
 *
 * This class may be used whenever you want to store or compute an alignment.
 * The underlying (ungapped sequence) remains unmodified, and is augmented by
 * gap information. Iterating over a gap_decorator sequences behaves as if you
 * had a rolled out, aligned sequence with alphabet and gap symbols.
 *
 * \details
 * The gap_decorator_anchor_set is a modified anchor list approach -- instead
 * of storing tuples of anchor positions and gap lengths relative to the
 * underlying sequence position, the anchor addresses are virtual and gap
 * lengths accumulated from left to right, i.e. an anchor gap stores at the 2nd
 * position all previous gap lengths including its own.
 * This reduces the lookup time to log (k) with k
 * being the number of continuous gaps (not gap symbols) and comes at the price
 * of O(k) worst-case runtime for gap modifications - when inserting or erasing
 * gaps the tailing gaps have to be updated by the resulting offset.
 * The gap maintaining structure is a set that is given a gap comparing structure.
 * Sets are implemented as red-black trees and perform random read operations in
 * log (k) time. The anchor set approach provides a good trade-off when using
 * both operator[] and gap insertion/erasure.
 *
 */
template <std::ranges::RandomAccessRange inner_type>
//!\cond
    requires container_concept<inner_type>
//!\endcond
class gap_decorator_anchor_set
{
    //!\brief The gap type as a tuple storing position and accumulated gap lengths.
    using gap_t = typename std::pair<size_t, size_t>;

    /*!\brief Structure allowing the comparison of gaps.
     * \tparam        gap_t Type of the gap storing position an length accumulator.
     * \details It is assumed that the gap structure is always in a consistent
     * state, i.e. there are no two gaps that are overlapping. The ordering of gaps
     * structures is exclusively dependent on the gap starting position which is
     * the ordering criterion for the anchor set implemented as ordered red-black tree.
     */
    template <typename gap_t>
    struct gap_compare {
        /*!\brief The operator allowing gap_t comparison.
         * \param[in] lhs   The left-hand gap to compare.
         * \param[in] rhs   The right-hand side gap to be compared.
         * \returns A boolean flag indicating whether the gap starting position
         *          of \p lhs is smaller than the one of \p rhs.
         */
         bool operator() (const gap_t& lhs, const gap_t& rhs) const {
             return lhs.first < rhs.first;
         }
     };
    //!\brief The iterator type for an anchor set.
    using set_iterator = typename std::set<gap_t, gap_compare<gap_t>>::iterator;

public:
    /*!\name Member types
     * \{
     */
    //!\brief The alphabet type of the underlying sequence.
    //!\hideinitializer
    using alphabet_type = typename ranges::v3::value_type_t<inner_type>;
    //!\brief The union type of the alphabet type and gap symbol type.
    //!\hideinitializer
    using value_type = gapped<alphabet_type>;
    //!\brief Use the value type as reference type.
    //!\hideinitializer
    using reference = value_type;
    //!\brief Use the const value type as reference type.
    //!\hideinitializer
    using const_reference = value_type const;
    //!\brief Use the size_type of the underlying sequence.
    //!\hideinitializer
    using size_type = size_type_t<inner_type>;
    //!\brief Use the difference_type of the underlying sequence.
    //!\hideinitializer
    using difference_type = typename ranges::v3::difference_type_t<inner_type>;
    //!\brief The iterator type of this container (a random access iterator).
    //!\hideinitializer
    using iterator = detail::random_access_iterator<gap_decorator_anchor_set>;
    //!\hideinitializer
    using const_iterator = detail::random_access_iterator<gap_decorator_anchor_set const>;
    //!\}

    /*!\name Constructors/Destructors
     * \{
     */
    //!\brief Default constructor.
    constexpr gap_decorator_anchor_set() = default;
    //!\brief Copy constructor.
    constexpr gap_decorator_anchor_set(gap_decorator_anchor_set const &) = default;
    //!\brief Copy construction via assignment.
    constexpr gap_decorator_anchor_set & operator=(gap_decorator_anchor_set const &) = default;
    //!\brief Move constructor.
    constexpr gap_decorator_anchor_set (gap_decorator_anchor_set && rhs) = default;
    //!\brief Move assignment.
    constexpr gap_decorator_anchor_set & operator=(gap_decorator_anchor_set && rhs) = default;
    //!\brief Construct by host sequence.
    constexpr gap_decorator_anchor_set(inner_type const & sequence): sequence{&sequence} {};
    //!\brief Direct sequence assignment resets previously inserted gaps.
    constexpr gap_decorator_anchor_set & operator=(inner_type const & sequence)
    {
        sequence = sequence;
        anchors.clear();
    };
    //!\brief Use default deconstructor.
    ~gap_decorator_anchor_set() = default;
    //!\}

    /*!\brief Returns the total length of the aligned sequence.
     * \returns The total length of the aligned sequence.
     *
     * \par Complexity
     *
     * Constant.
     *
     * \par Exceptions
     *
     * No-throw guarantee.
     */
    size_type size() const noexcept
    {
        if (anchors.size())
            return (*(anchors.rbegin())).second + sequence->size();
        return sequence->size();
    }

    /*!\brief Insert a gap of length size at the aligned sequence iterator position.
     * \param it    Iterator indicating the gap start position in the aligned sequence.
     * \param size  Number of gap symbols to be inserted.
     * \returns     An iterator pointing to the start position of the insertion.
     *
     * \par Complexity
     *
     * Average and worst case (insertion before last gap): o(k),
     * Best case (back insertion): o(log k).
     *
     * \par Exceptions
     * Throws assertion if \p it points beyond the end position.
     */
    iterator insert_gap(iterator const it, size_type const size = 1)
    {
        if (!size)
            return it;
        size_type const pos = it - begin();
        assert(pos <= this->size());

        set_iterator it_set = anchors.begin();
        // case 1: extend previous/surrounding gap already existing
        if ((pos < this->size()) && (((value_type)(*this)[pos] == gap::GAP) || (pos > 0 && (*this)[pos-1] == gap::GAP)))
        {
            it_set = anchors.lower_bound(gap_t{pos, 0/*Unused*/});
            if (it_set == anchors.end() || (*it_set).first > pos)
                it_set = std::prev(it_set);
            gap_t gap{(*it_set).first, (*it_set).second + size};
            // merge with successor
            auto it_next = it_set;
            ++it_next;
            if ((*it_set) < (*std::prev(anchors.end())) && (*it_next).first <= (*it_set).first + size - 1)
            {
                // extend gap for *it_next, delete *(it_next+1)
                gap.second += (*it_next).second;
                anchors.erase(it_next);
            }
            anchors.erase(it_set);
            anchors.insert(gap);
        }
        // case 2: create new anchor gap
        else
        {
            gap_t gap{pos, size};
            // pre: pos not in anchor set, find preceeding gap to add accumulated gaps
            if (anchors.size())
            {
                auto it_aux = anchors.lower_bound(gap_t{pos, 0/*Unused*/});
                if (it_aux != anchors.begin())
                    gap.second += (*--it_aux).second;
            }
            anchors.insert(gap);
        }
        // post-processing: reverse update of succeeding gaps
        rupdate(pos, size);
        return it;
    }

   /*!\brief Erase one gap symbol at the indicated iterator postion.
    * \param it     Iterator indicating the gap to be erased.
    * \returns      An iterator pointing to starting position of the gap erasure.
    *
    * \par Complexity
    *
    * O(log k)
    *
    * \par Exceptions
    * \throws seqan3::gap_erase_failure if character is no seqan3::gap.
    */
    iterator erase_gap(iterator const it)
    {
        // check if [it, it+gap_len[ covers [first, last[
        if ((*it) != gap::GAP) // [[unlikely]]
            throw gap_erase_failure("The range to be erased does not corresponds to a consecutive gap.");
        return erase_gap(it, it + 1);
    }

    /*!\brief Erase gap symbols at the iterator postions [first, last[.
     * \param[in]   first    The iterator pointing to the position where to start inserting gaps.
     * \param[in]   last     The iterator pointing to the position where to stop erasing gaps.
     * \returns     An iterator pointing to starting position of the gap erasure.
     *
     * \par Complexity
     *
     * O(log k)
     *
     * \par Exceptions
     *
     * \throws seqan3::gap_erase_failure if [\p first, \p last[ does not correspond
     * to a consecutive seqan3::gap range.
     */
    iterator erase_gap(iterator const first, iterator const last)
    {
        size_type pos1 = first - begin(), pos2 = last - begin();
        set_iterator it = anchors.lower_bound(gap_t{pos1, 0/*Unused*/});
        size_type gap_len = get_gap_length(it);

        if (it == anchors.end() || (*it).first > pos1)
            it = std::prev(it);
        // check if [it, it+gap_len[ covers [first, last[
        if (!(((*it).first <= pos1) && (((*it).first + gap_len) >= pos2))) // [[unlikely]]
            throw gap_erase_failure("The range to be erased does not corresponds to a consecutive gap.");
        // case 1: complete gap is deleted
        if (((*it).first == pos1) && (gap_len == pos2-pos1))
            anchors.erase(it);
        // case 2: gap to be deleted in tail or larger than 1 (equiv. to shift tail left, i.e. pos remains unchanged)
        else
        {
            gap_t gap{(*it).first, (*it).second - pos2 + pos1};
            anchors.erase(it);
            anchors.insert(gap);
        }
        // post-processing: forward update of succeeding gaps
        update(pos1, pos2-pos1);
        return first;
    }

    /*!\name Iterators
     * \{
     */
    /*!\brief Returns an iterator to the first element of the container.
     * \returns Iterator to the first element.
     *
     * If the container is empty, the returned iterator will be equal to end().
     *
     * \par Complexity
     *
     * Constant.
     *
     * \par Exceptions
     *
     * No-throw guarantee.
     */
    iterator begin() noexcept
    {
        return iterator{*this};
    }

    //!\copydoc begin()
    const_iterator begin() const noexcept
    {
        return const_iterator{*this};
    }

    //!\copydoc begin()
    const_iterator cbegin() const noexcept
    {
        return const_iterator{*this};
    }

    /*!\brief Returns an iterator to the element following the last element of the decorator.
     * \returns Iterator to the behind last element.
     *
     * This element acts as a placeholder; attempting to dereference it results in undefined behaviour.
     *
     * \par Complexity
     *
     * Constant.
     *
     * \par Exceptions
     *
     * No-throw guarantee.
     */
    iterator end() noexcept
    {
        return iterator{*this, size()};
    }

    //!\copydoc end()
    const_iterator end() const noexcept
    {
        return const_iterator{*this, size()};
    }

    //!\copydoc end()
    const_iterator cend() const noexcept
    {
        return const_iterator{*this, size()};
    }
    //!\}

    /*!\name Element access
     * \{
     */
    /*!\brief Return the i-th element as a reference.
     * \param i     The element to retrieve.
     * \returns     A reference of the gapped alphabet type.
     *
     * \par Complexity
     *
     * Constant.
     *
     * \par Exceptions
     *
     * Throws std::out_of_range exception if \p i is out of range.
     */
    reference at(size_type const i)
    {
        //TODO add SEQAN_UNLIKELY
        if (i >= size())
            throw std::out_of_range{"Trying to access element behind the last in gap_decorator."};
        return (*this)[i];
    }

    //!\copydoc at()
    const_reference at(size_type const i) const
    {
        //TODO add SEQAN_UNLIKELY
        if (i >= size())
            throw std::out_of_range{"Trying to access element behind the last in gap_decorator."};
        return (*this)[i];
    }

    //!\copydoc at()
    constexpr reference operator[](size_type const i) const
    {
        assert(i < size());
        // case 1: there are no gaps
        if (!anchors.size())
            return value_type((*sequence)[i]);
        // case 2: there are gaps
        set_iterator it = anchors.upper_bound(gap_t{i, 0/*Unused*/});
        if (it == anchors.begin())
            return value_type((*sequence)[i]); // since no gaps happen before i
        it = std::prev(it);
        size_type gap_len{(*it).second};
        if (it != anchors.begin())
            gap_len -= (*(std::prev(it, 1))).second;
        if (i < (*it).first + gap_len)
           return gap::GAP;
        else
           return value_type((*sequence)[i - (*it).second]);
    }
    //!\}

    /*!\name Comparison operators
     * \{
     */
    /*!\brief Compare gap decorators by underlying sequence and gaps.
     * \param[in] rhs   The right-hand side gap decorator to compare.
     * \returns A boolean flag indicating (in)equality of the aligned sequences.
     *
     * |par Complexity
     * Worst case: O(n*log k)
     *
     * \par Exceptions
     *
     * No-throw guarantee. Does not modify the aligned sequences.
     */
    template <typename inner_type2>
    //!\cond
        requires std::is_same_v<std::remove_const_t<inner_type>, std::remove_const_t<inner_type2>>
    //!\endcond
    constexpr bool operator==(gap_decorator_anchor_set<inner_type2> const & rhs) const noexcept
    {
        if (anchors.size() != rhs.anchors.size() || sequence->size() != rhs.sequence->size())
            return false;

        return std::ranges::equal(*this, rhs);
    }

    //!\copydoc operator==
    template <typename inner_type2>
    //!\cond
        requires std::is_same_v<std::remove_const_t<inner_type>, std::remove_const_t<inner_type2>>
    //!\endcond
    constexpr bool operator!=(gap_decorator_anchor_set<inner_type2> const & rhs) const noexcept
    {
        return !(operator==(rhs));
    }
    //!\}
private:
    /*!\brief Helper function to compute the length of the gap indicated by the
     * input iterator.
     * \param[in] it    Iterator over the internal gap set.
     * \returns The gap length corresponding to the gap pointed at by \p it.
     * \details The length of a gap is difference of the accumulator pointed at
     * and the one of its predecessor (if existing).
     *
     * \par Exceptions
     * No-throw guarantee.
     */
    constexpr size_type get_gap_length(set_iterator it) const noexcept
    {
        if (it == anchors.begin())
            return (*it).second;
        return (*it).second - (*std::prev(it)).second;
    }

    /*!\brief Update all anchor gaps after the indicated position by adding an offset.
     * \param[in] pos   Gap index after which to perform the update.
     * \param[in] size  Offset to be added to the virtual gap positions and its accumulators.
     *
     * |par Complexity
     * Linear in the number of gaps.
     *
     * \details For not invalidating the iterator over the ordered set, the
     * update is done in reverse manner excluding the indicated gap.
     */
    void rupdate(size_type const pos, size_type const size)
    {
        size_type new_key, new_val;
        for (auto it = std::prev(anchors.end(), 1); (*it).first > pos;)
        {
            new_key = (*it).first + size;
            new_val = (*it).second + size;
            anchors.emplace_hint(it, gap_t{new_key, new_val});
            anchors.erase(*it--);
        }
    }

    /*!\brief Update all anchor gaps after indicated position by substracting an offset.
     * \param[in] pos   Gap index after which to perform the update.
     * \param[in] size  Offset to be removed from the virtual gap positions and its accumulators.
     *
     * |par Complexity
     * Linear in the number of gaps.
     *
     * \details For not invalidating the iterator over the ordered set, the
     * decreasing is done in a forward manner excluding the indicated gap.
     *
     * \par Exceptions
     * Throws assert when initial update position is out of range.
     */
    void update(size_type const pos, size_type const size)
    {
        assert(pos < size);
        auto it = anchors.lower_bound(gap_t{pos + size + 1, 0/*Unused*/});
        while (it != anchors.end())
        {
            gap_t gap{(*it).first - size, (*it).second - size};
            anchors.insert(gap);
            set_iterator it_next = std::next(it);
            anchors.erase(it);
            it = it_next;
        }
    }

    //!\brief Stores a pointer to the ungapped, underlying sequence.
    inner_type const * sequence{};
    //!\brief Set storing the anchor gaps.
    std::set<gap_t, gap_compare<gap_t>> anchors{};
};

} // namespace seqan3
