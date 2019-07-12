// -----------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2019, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2019, Knut Reinert & MPI für molekulare Genetik
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at: https://github.com/seqan/seqan3/blob/master/LICENSE.md
// -----------------------------------------------------------------------------------------------------

/*!\file
 * \brief Provides seqan3::gap_decorator.
 * \author Marie Hoffmann <marie.hoffmann AT fu-berlin.de>
 * \author Svenja Mehringer <svenja.mehringer AT fu-berlin.de>
 */

#pragma once

#include <limits>
#include <set>
#include <tuple>
#include <type_traits>

#include <seqan3/alignment/exception.hpp>
#include <seqan3/alphabet/concept.hpp>
#include <seqan3/alphabet/gap/gap.hpp>
#include <seqan3/alphabet/gap/gapped.hpp>
#include <seqan3/range/container/concept.hpp>
#include <seqan3/range/detail/random_access_iterator.hpp>
#include <seqan3/range/view/view_all.hpp>
#include <seqan3/std/algorithm>
#include <seqan3/std/ranges>

namespace seqan3
{

/*!\brief A gap decorator allows the annotation of sequences with gap symbols
 *        while leaving the underlying sequence unmodified.
 * \tparam inner_type The type of range that will be decorated with gaps; must model std::ranges::RandomAccessRange
 *                    and std::ranges::SizedRange.
 * \implements seqan3::AlignedSequence
 * \ingroup decorator
 *
 * \details
 *
 * This class may be used whenever you want to store or compute an alignment. The underlying (ungapped) sequence
 * remains unmodified, and is augmented with gap information. The seqan3::gap_decorator behaves just like a
 * vector over a gapped alphabet when iterating over it, inserting/erasing gaps or accessing a position. The only
 * difference lies in the performance and size overhead (see below).
 *
 * ### Performance
 *
 * **n** The length of the underlying sequence.
 * **k** The number of contiguous gaps (not gap symbols).
 * **l** The total number of gap symbols.
 *
 * |            | access next | random access    | gap insert/erase at end | gap insert/erase random | size overhead |
 * |------------|-------------|----------------- |-------------------------|-------------------------|---------------|
 * | decorator  | \f$O(1)\f$  | \f$O(\log(k))\f$ | \f$O(\log(k))\f$        | \f$O(k)\f$              | \f$O(k)\f$    |
 * | vector     | \f$O(1)\f$  | \f$O(1)\f$       | \f$O(1)\f$              | \f$O(n)\f$              | \f$O(n)\f$    |
 *
 * The *size overhead* refers to the space that is needed when using each of the data structures in addition to an
 * already existing ungapped sequence.
 *
 * ### Implementation details
 *
 * This decorator stores a std::set over tuples of `(pos, cumulative_size)` where every entry represents one
 * contiguous stretch of gaps. `pos` is the (virtual) insert position in the underlying range and `cumulative_size`
 * is the length of that contiguous stretch of gaps plus the length of all preceding elements.
 * Resolving random access requires logarithmic access into the set and inserting or removing a gap symbol additionally
 * entails updating all subsequent elements in the set to preserve correct cumulative sizes.
 *
 * ### The seqan3::gap_decorator::iterator type
 *
 * \attention The iterator of the seqan3::gap_decorator does not model the
 *            [Cpp17InputIterator](https://en.cppreference.com/w/cpp/named_req/InputIterator) requirements of the
 *            STL because dereferencing the iterator returns a proxy and no operator-> is provided.
 *            Note that it does model the std::ranges::InputIterator.
 *
 */
template <std::ranges::ViewableRange inner_type>
//!\cond
    requires std::ranges::RandomAccessRange<inner_type> && std::ranges::SizedRange<inner_type> &&
             (std::is_const_v<std::remove_reference_t<inner_type>> || std::ranges::View<inner_type>)
//!\endcond
class gap_decorator
{
private:
    /*!\brief The iterator that moves over the seqan3::gap_decorator.
     *
     * \details
     *
     * This iterator returns values when dereferenced, not references, i.e. it does not satisfy the semantic
     * requirements of [LegacyForwardIterator](https://en.cppreference.com/w/cpp/named_req/ForwardIterator). It does
     * model the C++20 std::BidirectionalIterator (and std::ForwardIterator implicitly).
     */
    class gap_decorator_iterator
    {
    private:
        //!\brief Pointer to the underlying container structure.
        typename std::add_pointer_t<gap_decorator const> host{nullptr};
        //!\brief Stores the virtual position index for the seqan3::gap_decorator.
        typename gap_decorator::size_type pos{0u};
        //!\brief Stores the physical position in the ungapped/underlying view.
        int64_t ungapped_view_pos{0}; // must be signed because we need this value to be -1 in case of leading gaps.
        //!\brief Stores the position (inkl. gaps) where the last (consecutive) gap that is still before the current
        //!       iterator position ends.
        typename gap_decorator::size_type left_gap_end{0};
        //!\brief A pointer to the current anchor gap node. Note that the current tuple value at position 0 is the
        //!       start of the right gap that is still behind the current iterator position.
        typename gap_decorator::set_iterator_type anchor_set_it{};
        //!\brief Caches whether the iterator points to a gap (true) or not (false).
        bool is_at_gap{true};

        //!\brief A helper function that performs the random access into the anchor set, updating all member variables.
        void jump(typename gap_decorator::size_type const new_pos)
        {
            assert(new_pos <= host->size());
            pos = new_pos;

            anchor_set_it = host->anchors.upper_bound(anchor_gap_t{pos, host->bound_dummy});
            ungapped_view_pos = pos;

            if (anchor_set_it != host->anchors.begin())
            {
                typename gap_decorator::set_iterator_type prev{std::prev(anchor_set_it)};
                size_type gap_len{prev->second};

                if (prev != host->anchors.begin())
                    gap_len -= std::prev(prev)->second;

                ungapped_view_pos -= prev->second;
                left_gap_end = prev->first + gap_len;
            }

            if (ungapped_view_pos != static_cast<int64_t>(host->ungapped_view.size()) &&
                pos >= left_gap_end && (anchor_set_it == host->anchors.end() || pos < anchor_set_it->first))
                is_at_gap = false;
            else
                is_at_gap = true;
        }

    public:
        /*!\name Member types
         * \brief Make the parent's member types visible.
         * \{
         */
        //!\brief Type for distances between iterators.
        using difference_type = typename gap_decorator::difference_type;
        //!\brief Value type of container elements.
        using value_type = typename gap_decorator::value_type;
        //!\brief Const-reference type defined by container (which equals the reference type).
        using reference = typename gap_decorator::const_reference; // = reference
        //!\brief Equals reference type.
        using const_reference = reference;
        //!\brief Pointer of container value type.
        using pointer = value_type *;
        //!\brief Tag this class as a bidirectional iterator since random access is possible but not in constant time.
        using iterator_category = std::bidirectional_iterator_tag;
        //!\}

        /*!\name Constructors/Destructors
         * \{
         */
        //!\brief Default constructor.
        constexpr gap_decorator_iterator() = default;
        //!\brief Copy constructor.
        constexpr gap_decorator_iterator(gap_decorator_iterator const &) = default;
        //!\brief Copy construction via assignment.
        constexpr gap_decorator_iterator & operator=(gap_decorator_iterator const &) = default;
        //!\brief Move constructor.
        constexpr gap_decorator_iterator (gap_decorator_iterator &&) = default;
        //!\brief Move assignment.
        constexpr gap_decorator_iterator & operator=(gap_decorator_iterator &&) = default;
        //!\brief Use default deconstructor.
        ~gap_decorator_iterator() = default;

        //!\brief Construct from seqan3::gap_decorator and initialise members.
        explicit constexpr gap_decorator_iterator(gap_decorator const & host_) noexcept :
            host(&host_), anchor_set_it{host_.anchors.begin()}
        {
            if (host_.anchors.size() && (*host_.anchors.begin()).first == 0) // there are gaps at the very front
            {
                --ungapped_view_pos; // set ungapped_view_pos to -1 so operator++ works without an extra if-branch.
                left_gap_end = anchor_set_it->second;
                ++anchor_set_it;
            }
            else
            {
                is_at_gap = false;
            }
        }

        //!\brief Construct from seqan3::gap_decorator and explicit position.
        constexpr gap_decorator_iterator(gap_decorator const & host_,
                                                    typename gap_decorator::size_type const pos_) noexcept :
             host(&host_)
        {
            jump(pos_); // random access to pos
        }
        //!\}

        /*!\name Arithmetic operators
         * \{
        */
        //!\brief Pre-increment, returns updated iterator.
        constexpr gap_decorator_iterator & operator++() noexcept
        {
            assert(host); // host is set
            ++pos;

            if (pos < left_gap_end)
            {   // we stay within the preceding gap stretch
                // no op but must precede other cases
            }
            else if (anchor_set_it == host->anchors.end() || pos < anchor_set_it->first)
            {   // proceed within the view since we are right of the previous gap but didn't arrive at the right gap yet
                ++ungapped_view_pos;
                if (ungapped_view_pos != static_cast<int64_t>(host->ungapped_view.size()))
                    is_at_gap = false;
            }
            else
            {   // we arrived at the right gap and have to update the variables. ungapped_view_pos remains unchanged.
                left_gap_end = anchor_set_it->first + anchor_set_it->second -
                               ((anchor_set_it != host->anchors.begin()) ? (std::prev(anchor_set_it))->second : 0);
                ++anchor_set_it;
                is_at_gap = true;
            }

            return *this;
        }

        //!\brief Pre-decrement, returns updated iterator.
        constexpr gap_decorator_iterator & operator--() noexcept
        {
            assert(host); // host is set
            --pos;

            if (pos < left_gap_end)
            {   // there was no gap before but we arrive at the left gap and have to update the variables.
                (anchor_set_it != host->anchors.begin()) ? --anchor_set_it : anchor_set_it;

                if (anchor_set_it != host->anchors.begin())
                {
                    auto prev = std::prev(anchor_set_it);
                    left_gap_end = prev->first + prev->second -
                                   ((prev != host->anchors.begin()) ? std::prev(prev)->second : 0);
                }
                else // [[unlikely]]
                {
                    left_gap_end = 0;
                }
                is_at_gap = true;
            }
            else if (anchor_set_it == host->anchors.end() || pos < anchor_set_it->first)
            {   // we are neither at the left nor right gap
                --ungapped_view_pos;
                is_at_gap = false;
            }
            // else -> no op (we are still within the right gap stretch)

            return *this;
        }

        //!\brief Post-increment, returns previous iterator state (delegates to pre-increment).
        constexpr gap_decorator_iterator operator++(int) noexcept
        {
            gap_decorator_iterator cpy{*this};
            ++(*this);
            return cpy;
        }

        //!\brief Post-decrement, returns previous iterator state (delegates to pre-decrement).
        constexpr gap_decorator_iterator operator--(int) noexcept
        {
            gap_decorator_iterator cpy{*this};
            --(*this);
            return cpy;
        }
        //!\}

        /*!\name Reference/Dereference operators
         * \{
        */
        //!\brief Dereference operator returns a copy of the element currently pointed at.
        constexpr reference operator*() const noexcept
        {
            return (is_at_gap) ? static_cast<reference>(gap{})
                               : static_cast<reference>(host->ungapped_view[ungapped_view_pos]);
        }
        //!\}

        /*!\name Comparison operators
         * \brief Compares iterators by virtual position.
         * \{
         */

        //!\brief Checks whether `*this` is equal to `rhs`.
        constexpr friend bool operator==(gap_decorator_iterator const & lhs,
                                         gap_decorator_iterator const & rhs)
        {
            return lhs.pos == rhs.pos;
        }

        //!\brief Checks whether `*this` is not equal to `rhs`.
        constexpr friend bool operator!=(gap_decorator_iterator const & lhs,
                                         gap_decorator_iterator const & rhs)
        {
            return lhs.pos != rhs.pos;
        }

        //!\brief Checks whether `*this` is less than `rhs`.
        constexpr friend bool operator<(gap_decorator_iterator const & lhs,
                                        gap_decorator_iterator const & rhs)
        {
            return lhs.pos < rhs.pos;
        }

        //!\brief Checks whether `*this` is greater than `rhs`.
        constexpr friend bool operator>(gap_decorator_iterator const & lhs,
                                        gap_decorator_iterator const & rhs)
        {
            return lhs.pos > rhs.pos;
        }

        //!\brief Checks whether `*this` is less than or equal to `rhs`.
        constexpr friend bool operator<=(gap_decorator_iterator const & lhs,
                                         gap_decorator_iterator const & rhs)
        {
            return lhs.pos <= rhs.pos;
        }

        //!\brief Checks whether `*this` is greater than or equal to `rhs`.
        constexpr friend bool operator>=(gap_decorator_iterator const & lhs,
                                         gap_decorator_iterator const & rhs)
        {
            return lhs.pos >= rhs.pos;
        }
        //!\}
    };

public:
    /*!\name Range-associated member types
     * \{
     */
    //!\brief The variant type of the alphabet type and gap symbol type (see seqan3::gapped).
    using value_type = gapped<value_type_t<inner_type>>;
    //!\brief Use the value type as reference type because the underlying sequence must not be modified.
    using reference = value_type;
    //!\brief const_reference type equals reference type equals value type because the underlying sequence must not
    //!       be modified.
    using const_reference = reference;
    //!\brief The size_type of the underlying sequence.
    using size_type = size_type_t<inner_type>;
    //!\brief The difference type of the underlying sequence.
    using difference_type = difference_type_t<inner_type>;
    //!\brief The iterator type of this container (a bidirectional iterator).
    using iterator = gap_decorator_iterator;
    //!\brief The const_iterator equals the iterator type. Since no references are ever returned and thus the underlying
    //!        sequence cannot be modified through the iterator there is no need for const.
    using const_iterator = iterator;
    //!\}

    //!\brief The underlying ungapped range type.
    using unaligned_seq_type = inner_type;

    /*!\name Constructors, destructor and assignment.
     * \{
     */
    //!\brief Default constructor. Attention: all operations on a solely default constructed decorator,
    //!       except assigning a new range, are UB.
    constexpr gap_decorator() = default;
    //!\brief Copy constructor.
    constexpr gap_decorator(gap_decorator const &) = default;
    //!\brief Copy construction via assignment.
    constexpr gap_decorator & operator=(gap_decorator const &) = default;
    //!\brief Move constructor.
    constexpr gap_decorator(gap_decorator && rhs) = default;
    //!\brief Move assignment.
    constexpr gap_decorator & operator=(gap_decorator && rhs) = default;
    //!\brief Use default deconstructor.
    ~gap_decorator() = default;

    //!\brief Construct with the ungapped range type.
    template <typename other_range_t>
    //!\cond
         requires !std::Same<other_range_t, gap_decorator> &&
                  std::Same<remove_cvref_t<other_range_t>, remove_cvref_t<inner_type>> &&
                  std::ranges::ViewableRange<other_range_t> // at end, otherwise it competes with the move ctor
    //!\endcond
    gap_decorator(other_range_t && range) : ungapped_view{view::all(std::forward<inner_type>(range))}
    {} // TODO (@smehringer) only works for copyable views. Has to be changed once views are not required to be copyable anymore.
    // !\}

    /*!\brief Returns the total length of the aligned sequence.
     * \returns The total length of the aligned sequence (gaps included).
     *
     * ### Complexity
     *
     * Constant.
     *
     * ### Exceptions
     *
     * No-throw guarantee.
     */
    size_type size() const noexcept
    {
        if (anchors.size())
            return anchors.rbegin()->second + ungapped_view.size();

        return ungapped_view.size();
    }

    /*!\name Aligned sequence modifications
     * \{
     */
    /*!\brief Insert a gap of length count at the aligned sequence iterator position.
     * \param it     Iterator indicating the gap start position in the aligned sequence.
     * \param count  Number of gap symbols to be inserted.
     * \returns      An iterator pointing to the start position of the insertion.
     *
     * ### Complexity
     *
     * Average and worst case (insertion before last gap): \f$O(k)\f$,
     * Best case (back insertion): \f$O(\log k)\f$.
     */
    iterator insert_gap(iterator const it, size_type const count = 1)
    {
        if (!count) // [[unlikely]]
            return it;

        size_type const pos = std::distance(begin(), it);
        assert(pos <= size());

        set_iterator_type it_set = anchors.upper_bound(anchor_gap_t{pos, bound_dummy});

        if (it_set == anchors.begin()) // will also catch if anchors is empty since begin() == end()
        {
            anchors.emplace_hint(anchors.begin(), anchor_gap_t{pos, count});
        }
        else // there are gaps before pos
        {
            --it_set;
            auto gap_len{it_set->second};
            if (it_set != anchors.begin())
                gap_len -= (*(std::prev(it_set))).second;

            if (it_set->first + gap_len >= pos) // extend existing gap
            {
                anchor_gap_t gap{it_set->first, it_set->second + count};
                it_set = anchors.erase(it_set);
                anchors.insert(it_set, gap);
            }
            else                                  // insert new gap
            {
                anchor_gap_t gap{pos, it_set->second + count};
                ++it_set;
                anchors.insert(it_set, gap);
            }
        }

        // post-processing: reverse update of succeeding gaps
        rupdate(pos, count);
        return iterator{*this, pos};
    }

   /*!\brief Erase one gap symbol at the indicated iterator postion.
    * \param it     Iterator indicating the gap to be erased.
    * \returns      Iterator following the last removed element.
    * \throws seqan3::gap_erase_failure if character is no seqan3::gap.
    *
    * \details
    *
    * ### Complexity
    *
    * \f$O(\log k)\f$
    */
    iterator erase_gap(iterator const it)
    {
        // check if [it, it+gap_len[ covers [first, last[
        if ((*it) != gap{}) // [[unlikely]]
            throw gap_erase_failure("The range to be erased does not correspond to a consecutive gap.");

        auto end_it = std::next(it);
        return erase_gap(it, end_it);
    }

    /*!\brief Erase gap symbols at the iterator postions [first, last[.
     * \param[in]   first    The iterator pointing to the position where to start inserting gaps.
     * \param[in]   last     The iterator pointing to the position where to stop erasing gaps.
     * \returns     Iterator following the last removed element.
     * \throws seqan3::gap_erase_failure if [\p first, \p last[ does not correspond
     * to a consecutive range of seqan3::gap 's.
     *
     * \details
     *
     * ### Complexity
     *
     * \f$O(\log k)\f$
     */
    iterator erase_gap(iterator const first, iterator const last)
    {
        size_type const pos1 = std::distance(begin(), first);
        size_type const pos2 = std::distance(begin(), last);
        set_iterator_type it = anchors.upper_bound(anchor_gap_t{pos1, bound_dummy}); // first element greater than pos1

        if (it == anchors.begin())
            throw gap_erase_failure{"There is no gap to erase in range [" + std::to_string(pos1) + "," +
                                    std::to_string(pos2) + "]."};

        --it;
        size_type const gap_len = gap_length(it);

        // check if [it, it+gap_len[ covers [first, last[
        if ((it->first + gap_len) < pos2) // [[unlikely]]
        {
            throw gap_erase_failure{"The range to be erased does not correspond to a consecutive gap."};
        }
        // case 1: complete gap is deleted
        else if (gap_len == pos2 - pos1)
        {
            it = anchors.erase(it);
        }
        // case 2: gap to be deleted in tail or larger than 1 (equiv. to shift tail left, i.e. pos remains unchanged)
        else
        {
            anchor_gap_t gap{it->first, it->second - pos2 + pos1};
            it = anchors.erase(it);
            it = anchors.insert(it, gap); // amortized constant because of hint
            ++it; // update node after the current
        }

        // post-processing: forward update of succeeding gaps
        update(it, pos2 - pos1);

        return iterator{*this, pos1};
    }

    /*!\brief Assigns a new sequence of type seqan3::gap_decorator::unaligned_seq_type to the decorator.
     * \param[in,out] dec       The decorator to modify.
     * \param[in]     unaligned The unaligned sequence to assign.
     */
    template <typename unaligned_seq_t> // generic template to use forwarding reference
    //!\cond
        requires std::Assignable<gap_decorator &, unaligned_seq_t>
    //!\endcond
    friend void assign_unaligned(gap_decorator & dec, unaligned_seq_t && unaligned)
    {
        dec = unaligned;
    }
    //!\}

    /*!\name Iterators
     * \{
     */
    /*!\brief Returns an iterator to the first element of the container.
     * \returns Iterator to the first element.
     *
     * If the container is empty, the returned iterator will be equal to end().
     *
     * ### Complexity
     *
     * Constant.
     *
     * ### Exceptions
     *
     * No-throw guarantee.
     */
    iterator begin() const noexcept
    {
        return iterator{*this};
    }

    //!\copydoc begin()
    const_iterator cbegin() const noexcept
    {
        return const_iterator{*this};
    }

    /*!\brief Returns an iterator to the element following the last element of the decorator.
     * \returns Iterator to the behind last element.
     *
     * \attention This element acts as a placeholder; attempting to dereference it results in undefined behaviour.
     *
     * ### Complexity
     *
     * Constant.
     *
     * ### Exceptions
     *
     * No-throw guarantee.
     */
    iterator end() const noexcept
    {
        return iterator{*this, size()};
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
     * ### Complexity
     *
     * \f$O(\log k)\f$ where \f$k\f$ is the number of gaps.
     *
     * ### Exceptions
     *
     * Throws std::out_of_range exception if \p i is out of range.
     */
    reference at(size_type const i)
    {
        if (i >= size()) // [[unlikely]]
            throw std::out_of_range{"Trying to access element behind the last in gap_decorator."};
        return (*this)[i];
    }

    //!\copydoc at()
    const_reference at(size_type const i) const
    {
        if (i >= size()) // [[unlikely]]
            throw std::out_of_range{"Trying to access element behind the last in gap_decorator."};
        return (*this)[i];
    }

    /*!\brief Return the i-th element as a reference.
     * \param i     The element to retrieve.
     * \returns     A reference of the gapped alphabet type.
     *
     * This function delegates to an iterator seqan3::gap_decorator.
     *
     * ### Complexity
     *
     * \f$O(\log k)\f$ where \f$k\f$ is the number of gaps.
     *
     */
    constexpr reference operator[](size_type const i) const noexcept
    {
        return *iterator{*this, i};
    }
    //!\}

    /*!\name Comparison operators
     * \brief Compares two seqan3::gap_decorator 's by underlying sequence and gaps.
     * \param[in] lhs The left-hand side gap decorator to compare.
     * \param[in] rhs The right-hand side gap decorator to compare.
     * \returns A boolean flag indicating (in)equality of the aligned sequences.
     *
     * \details
     *
     * ### Complexity
     *
     * Worst case: \f$O(n*\log k)\f$
     * Constant in case the decorators have not the same number of (consecutive) gaps.
     *
     * ### Exceptions
     *
     * No-throw guarantee. Does not modify the aligned sequences.
     * \{
     */

    //!\brief Checks whether `lhs` is equal to `rhs`.
    friend bool operator==(gap_decorator const & lhs, gap_decorator const & rhs) noexcept
    {
        if (lhs.size()  == rhs.size()  &&
            lhs.anchors == rhs.anchors &&
            std::ranges::equal(lhs.ungapped_view, rhs.ungapped_view))
        {
            return true;
        }

        return false;
    }

    //!\brief Checks whether `lhs` is not equal to `rhs`.
    friend bool operator!=(gap_decorator const & lhs, gap_decorator const & rhs) noexcept
    {
        return !(lhs == rhs);
    }

    //!\brief Checks whether `lhs` is less than `rhs`.
    friend bool operator<(gap_decorator const & lhs, gap_decorator const & rhs) noexcept
    {
        auto lit = lhs.begin();
        auto rit = rhs.begin();

        while (lit != lhs.end() && rit != rhs.end() && *lit == *rit)
            ++lit, ++rit;

        if (rit == rhs.end())
            return false;           //  lhs == rhs, or rhs prefix of lhs
        else if (lit == lhs.end())
            return true;            // lhs prefix of rhs

        return *lit < *rit;
    }

    //!\brief Checks whether `lhs` is less than or equal to `rhs`.
    friend bool operator<=(gap_decorator const & lhs, gap_decorator const & rhs) noexcept
    {
        auto lit = lhs.begin();
        auto rit = rhs.begin();

        while (lit != lhs.end() && rit != rhs.end() && *lit == *rit)
            ++lit, ++rit;

        if (lit == lhs.end())
            return true;            // lhs == rhs, or lhs prefix of rhs
        else if (rit == rhs.end())
            return false;           // rhs prefix of lhs

        return *lit < *rit;
    }

    //!\brief Checks whether `lhs` is greater than `rhs`.
    friend bool operator>(gap_decorator const & lhs, gap_decorator const & rhs) noexcept
    {
        return !(lhs <= rhs);
    }

    //!\brief Checks whether `lhs` is greater than or equal to `rhs`.
    friend bool operator>=(gap_decorator const & lhs, gap_decorator const & rhs) noexcept
    {
        return !(lhs < rhs);
    }
    //!\}

private:
    //!\brief The gap type as a tuple storing position and accumulated gap lengths.
    using anchor_gap_t = typename std::pair<size_t, size_t>;

    //!\brief The type of set to store the anchor gaps.
    using anchor_set_type = std::set<anchor_gap_t>;

    //!\brief The iterator type for an anchor set.
    using set_iterator_type = typename anchor_set_type::iterator;

    //!\brief The maximum value is needed for a correct search with upper_bound() in the anchor set.
    constexpr static size_t bound_dummy{std::numeric_limits<size_t>::max()};

    /*!\brief Helper function to compute the length of the gap indicated by the input iterator.
     * \param[in] it    Iterator over the internal gap set.
     * \returns The gap length corresponding to the gap pointed at by \p it.
     *
     * \details
     *
     * The length of a gap pointed to by \p it is the difference of the current cumulative sum of gaps
     * (second tuple position) and the one of its predecessor (if existing).
     *
     * ### Exceptions
     * No-throw guarantee.
     */
    constexpr size_type gap_length(set_iterator_type it) const noexcept
    {
        return (it == anchors.begin()) ? it->second : it->second - (*std::prev(it)).second;
    }

    /*!\brief Update all anchor gaps after the indicated position by adding an offset.
     * \param[in] pos    Gap index after which to perform the update.
     * \param[in] offset Offset to be added to the virtual gap positions and its accumulators.
     *
     * \details
     *
     * In order to avoid key conflicts when inserting into anchors, the
     * update is done in reverse manner excluding the indicated gap.
     *
     * ### Complexity
     * Linear in the number of gaps.
     */
    void rupdate(size_type const pos, size_type const offset)
    {
        for (auto it = std::prev(anchors.end(), 1); it->first > pos;)
        {
            anchors.emplace_hint(it, anchor_gap_t{it->first + offset, it->second + offset});
            anchors.erase(*it--);
        }
    }

    /*!\brief Update all anchor gaps after indicated position by substracting an offset.
     * \param[in] it     Iterator pointing to the first anchor node to update.
     * \param[in] offset Offset to be removed from the virtual gap positions and its accumulators.
     *
     * \details
     *
     * In order to avoid key conflicts when inserting into anchors, the
     * decreasing is done in a forward manner excluding the indicated gap.
     *
     * ### Complexity
     * Linear in the number of gaps.
     */
    void update(set_iterator_type it, size_type const offset)
    {
        while (it != anchors.end())
        {
            anchor_gap_t gap{it->first - offset, it->second - offset};
            it = anchors.erase(it);
            it = anchors.insert(it, gap);
            ++it;
        }
    }

    //!\brief Stores a (copy of a) view to the ungapped, underlying sequence.
    decltype(view::all(std::declval<inner_type &&>())) ungapped_view{};

    //!\brief Set storing the anchor gaps.
    anchor_set_type anchors{};
};

/*!\name Type deduction guides
 * \{
 */
//!\brief Ranges (not views!) always deduce to `const & range_type` since they are access-only anyway.
template <std::ranges::ViewableRange urng_t>
//!\cond
    requires !std::ranges::View<std::remove_reference_t<urng_t>>
//!\endcond
gap_decorator(urng_t && range) -> gap_decorator<std::remove_reference_t<urng_t> const &>;

//!\brief Views always deduce to their respective type because they are copied.
template <std::ranges::View urng_t>
gap_decorator(urng_t range) -> gap_decorator<urng_t>;
//!\}

} // namespace seqan

namespace seqan3::detail
{

//!\brief Type trait that declares any seqan3::gap_decorator to be **NOT a view**.
template <typename type>
constexpr int enable_view<seqan3::gap_decorator<type>> = 0;

template <typename type>
constexpr int enable_view<seqan3::gap_decorator<type> const> = 0;

} // namespace seqan3::detail
