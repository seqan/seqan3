// SPDX-FileCopyrightText: 2006-2025 Knut Reinert & Freie Universität Berlin
// SPDX-FileCopyrightText: 2016-2025 Knut Reinert & MPI für molekulare Genetik
// SPDX-License-Identifier: BSD-3-Clause

/*!\file
 * \brief Provides seqan3::gap_decorator.
 * \author Marie Hoffmann <marie.hoffmann AT fu-berlin.de>
 * \author Svenja Mehringer <svenja.mehringer AT fu-berlin.de>
 */

#pragma once

#include <algorithm>
#include <limits>
#include <ranges>
#include <set>
#include <tuple>
#include <type_traits>

#include <seqan3/alignment/exception.hpp>
#include <seqan3/alphabet/concept.hpp>
#include <seqan3/alphabet/gap/gap.hpp>
#include <seqan3/alphabet/gap/gapped.hpp>
#include <seqan3/utility/views/type_reduce.hpp>

namespace seqan3
{

/*!\brief A gap decorator allows the annotation of sequences with gap symbols
 *        while leaving the underlying sequence unmodified.
 * \tparam inner_type The type of range that will be decorated with gaps; must model std::ranges::random_access_range
 *                    and std::ranges::sized_range.
 * \implements seqan3::writable_aligned_sequence
 * \ingroup alignment_decorator
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
 *            [Cpp17BidirectionalIterator](https://en.cppreference.com/w/cpp/iterator/bidirectional_iterator)
 *            requirements of the STL because dereferencing the iterator returns a proxy and no operator-> is provided.
 *            It does model the C++20 std::bidirectional_iterator.
 *
 * \stableapi{Since version 3.1.}
 */
template <std::ranges::viewable_range inner_type>
    requires std::ranges::random_access_range<inner_type> && std::ranges::sized_range<inner_type>
          && (std::is_const_v<std::remove_reference_t<inner_type>> || std::ranges::view<inner_type>)
class gap_decorator
{
private:
    // Declaration of class's iterator types; for the definition see below.
    template <bool = true>
    class basic_iterator;

    //!\brief The iterator type of this container (a bidirectional iterator).
    using iterator = basic_iterator<true>;
    //!\brief The const_iterator equals the iterator type. Since no references are ever returned and thus the underlying
    //!        sequence cannot be modified through the iterator there is no need for const.
    using const_iterator = iterator;

    //!\brief The type of the underlying view wrapped in seqan3::views::type_reduce.
    using ungapped_view_type = decltype(views::type_reduce(std::declval<inner_type &&>()));

public:
    /*!\name Range-associated member types
     * \{
     */
    /*!\brief The variant type of the alphabet type and gap symbol type (see seqan3::gapped).
     * \details
     * \stableapi{Since version 3.1.}
     */
    using value_type = gapped<std::ranges::range_value_t<inner_type>>;

    /*!\brief Use the value type as reference type because the underlying sequence must not be modified.
     * \details
     * \stableapi{Since version 3.1.}
     */
    using reference = value_type;

    /*!\brief const_reference type equals reference type equals value type because the underlying sequence must not
     *        be modified.
     * \details
     * \stableapi{Since version 3.1.}
     */
    using const_reference = reference;

    /*!\brief The size_type of the underlying sequence.
     * \details
     * \stableapi{Since version 3.1.}
     */
    using size_type = std::ranges::range_size_t<inner_type>;

    /*!\brief The difference type of the underlying sequence.
     * \details
     * \stableapi{Since version 3.1.}
     */
    using difference_type = std::ranges::range_difference_t<inner_type>;
    //!\}

    /*!\brief The underlying ungapped range type.
     *
     * \stableapi{Since version 3.1.}
     */
    using unaligned_sequence_type = inner_type;

    /*!\name Constructors, destructor and assignment
     * \{
     */
    /*!\brief Default constructor.
     * \details
     * \attention All operations on a default constructed decorator, except assigning a new range
     *            (e.g., seqan3::gap_decorator::assign_unaligned), are undefined behaviour.
     *
     * \stableapi{Since version 3.1.}
     */
    gap_decorator() = default;
    gap_decorator(gap_decorator const &) = default;             //!< Defaulted.
    gap_decorator & operator=(gap_decorator const &) = default; //!< Defaulted.
    gap_decorator(gap_decorator && rhs) = default;              //!< Defaulted.
    gap_decorator & operator=(gap_decorator && rhs) = default;  //!< Defaulted.
    ~gap_decorator() = default;                                 //!< Defaulted.

    /*!\brief Construct with the ungapped range type.
     * \details
     * \experimentalapi{Experimental since version 3.1. This is a non-standard C++ extension.}
     */
    template <typename other_range_t>
        requires (!std::same_as<other_range_t, gap_decorator>)
              && std::same_as<std::remove_cvref_t<other_range_t>, std::remove_cvref_t<inner_type>>
              && std::ranges::viewable_range<other_range_t> // at end, otherwise it competes with the move ctor
    gap_decorator(other_range_t && range) : ungapped_view{views::type_reduce(std::forward<inner_type>(range))}
    {
    } // TODO (@smehringer) only works for copyable views. Has to be changed once views are not required to be copyable anymore.
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
     * Strong exception guarantee.
     *
     * \stableapi{Since version 3.1.}
     */
    size_type size() const
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
     *
     * \stableapi{Since version 3.1.}
     */
    iterator insert_gap(const_iterator const it, size_type const count = 1)
    {
        if (!count) // [[unlikely]]
            return it;

        size_type const pos = it - begin();
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
            else // insert new gap
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
     *
     * \stableapi{Since version 3.1.}
    */
    iterator erase_gap(const_iterator const it)
    {
        // check if [it, it+gap_len[ covers [first, last[
        if ((*it) != gap{}) // [[unlikely]]
            throw gap_erase_failure("The range to be erased does not correspond to a consecutive gap.");

        return erase_gap(it, std::next(it));
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
     *
     * \stableapi{Since version 3.1.}
     */
    iterator erase_gap(const_iterator const first, const_iterator const last)
    {
        size_type const pos1 = first - begin();
        size_type const pos2 = last - begin();
        set_iterator_type it = anchors.upper_bound(anchor_gap_t{pos1, bound_dummy}); // first element greater than pos1

        if (it == anchors.begin())
            throw gap_erase_failure{"There is no gap to erase in range [" + std::to_string(pos1) + ","
                                    + std::to_string(pos2) + "]."};

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
            ++it;                         // update node after the current
        }

        // post-processing: forward update of succeeding gaps
        update(it, pos2 - pos1);

        return iterator{*this, pos1};
    }

    /*!\brief Assigns a new sequence of type seqan3::gap_decorator::unaligned_sequence_type to the decorator.
     * \param[in,out] dec       The decorator to modify.
     * \param[in]     unaligned The unaligned sequence to assign.
     * \details
     *
     * \experimentalapi{Experimental since version 3.1.}
     */
    template <typename unaligned_sequence_t> // generic template to use forwarding reference
        requires std::assignable_from<gap_decorator &, unaligned_sequence_t>
    friend void assign_unaligned(gap_decorator & dec, unaligned_sequence_t && unaligned)
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
     *
     * \stableapi{Since version 3.1.}
     */
    const_iterator begin() const noexcept
    {
        return iterator{*this};
    }

    //!\copydoc begin()
    const_iterator cbegin() const noexcept
    {
        return const_iterator{*this};
    }

    /*!\brief Returns an iterator pointing behind the last element of the decorator.
     * \returns Iterator pointing behind the last element.
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
     *
     * \stableapi{Since version 3.1.}
     */
    const_iterator end() const noexcept
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
     *
     * \stableapi{Since version 3.1.}
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
     * \stableapi{Since version 3.1.}
     */
    reference operator[](size_type const i) const
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
     * Strong exception guarantee.
     * \{
     */

    /*!\brief Checks whether `lhs` is equal to `rhs`.
     * \details
     * \stableapi{Since version 3.1.}
     */
    friend bool operator==(gap_decorator const & lhs, gap_decorator const & rhs)
    {
        if (lhs.size() == rhs.size() && lhs.anchors == rhs.anchors
            && std::ranges::equal(lhs.ungapped_view, rhs.ungapped_view))
        {
            return true;
        }

        return false;
    }

    /*!\brief Checks whether `lhs` is not equal to `rhs`.
     * \details
     * \stableapi{Since version 3.1.}
     */
    friend bool operator!=(gap_decorator const & lhs, gap_decorator const & rhs)
    {
        return !(lhs == rhs);
    }

    /*!\brief Checks whether `lhs` is less than `rhs`.
     * \details
     * \stableapi{Since version 3.1.}
     */
    friend bool operator<(gap_decorator const & lhs, gap_decorator const & rhs)
    {
        auto lit = lhs.begin();
        auto rit = rhs.begin();

        while (lit != lhs.end() && rit != rhs.end() && *lit == *rit)
            ++lit, ++rit;

        if (rit == rhs.end())
            return false; //  lhs == rhs, or rhs prefix of lhs
        else if (lit == lhs.end())
            return true; // lhs prefix of rhs

        return *lit < *rit;
    }

    /*!\brief Checks whether `lhs` is less than or equal to `rhs`.
     * \details
     * \stableapi{Since version 3.1.}
     */
    friend bool operator<=(gap_decorator const & lhs, gap_decorator const & rhs)
    {
        auto lit = lhs.begin();
        auto rit = rhs.begin();

        while (lit != lhs.end() && rit != rhs.end() && *lit == *rit)
            ++lit, ++rit;

        if (lit == lhs.end())
            return true; // lhs == rhs, or lhs prefix of rhs
        else if (rit == rhs.end())
            return false; // rhs prefix of lhs

        return *lit < *rit;
    }

    /*!\brief Checks whether `lhs` is greater than `rhs`.
     * \details
     * \stableapi{Since version 3.1.}
     */
    friend bool operator>(gap_decorator const & lhs, gap_decorator const & rhs)
    {
        return !(lhs <= rhs);
    }

    /*!\brief Checks whether `lhs` is greater than or equal to `rhs`.
     * \details
     * \stableapi{Since version 3.1.}
     */
    friend bool operator>=(gap_decorator const & lhs, gap_decorator const & rhs)
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
    static constexpr size_t bound_dummy{std::numeric_limits<size_t>::max()};

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
     * Strong exception guarantee.
     */
    size_type gap_length(set_iterator_type it) const
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
    ungapped_view_type ungapped_view{};

    //!\brief Set storing the anchor gaps.
    anchor_set_type anchors{};
};

/*!\name Type deduction guides
 * \{
 */
/*!\brief Ranges (not views!) always deduce to `const & range_type` since they are access-only anyway.
 * \details
 * \experimentalapi{Experimental since version 3.1.}
 */
template <std::ranges::viewable_range urng_t>
    requires (!std::ranges::view<std::remove_reference_t<urng_t>>)
gap_decorator(urng_t && range) -> gap_decorator<std::remove_reference_t<urng_t> const &>;

/*!\brief Views always deduce to their respective type because they are copied.
 * \details
 * \experimentalapi{Experimental since version 3.1.}
 */
template <std::ranges::view urng_t>
gap_decorator(urng_t range) -> gap_decorator<urng_t>;
//!\}

/*!\brief The iterator type over a seqan3::gap_decorator.
 * \implements seqan3::pseudo_random_access_iterator
 *
 * \details
 *
 * This iterator returns values when dereferenced, not references, i.e. it does not satisfy the semantic
 * requirements of [Cpp17BidirectionalIterator](https://en.cppreference.com/w/cpp/iterator/bidirectional_iterator).
 * It does model the C++20 std::bidirectional_iterator. In addition, it offers all interfaces of a standard
 * std::random_access_iterator except the iterator category which is defined as std::bidirectional_iterator_tag,
 * because the complexity of the iterator is logarithmic and not constant. However, all interfaces inside the
 * seqan3::gap_decorator make use of the more efficient logarithmic implementation. Be aware, that if you want to use
 * the seqan3::gap_decorator in a generic algorithm, e.g. std::ranges::distance, the slower linear version will be
 * picked due to the constraints of the iterator category. To achieve optimal performance in a generic context you can
 * use seqan3::views::enforce_random_access to get a range over seqan3::gap_decorator which models
 * std::random_access_iterator albeit its non-conforming runtime complexity.
 */
template <std::ranges::viewable_range inner_type>
    requires std::ranges::random_access_range<inner_type> && std::ranges::sized_range<inner_type>
          && (std::is_const_v<std::remove_reference_t<inner_type>> || std::ranges::view<inner_type>)
template <bool>
class gap_decorator<inner_type>::basic_iterator
{
protected:
    //!\brief Pointer to the underlying container structure.
    typename std::add_pointer_t<gap_decorator const> host{nullptr};
    //!\brief Stores the virtual position index for the seqan3::gap_decorator.
    typename gap_decorator::size_type pos{0u};
    //!\brief Stores the physical position in the ungapped/underlying view.
    int64_t ungapped_view_pos{0}; // must be signed because we need this value to be -1 in case of leading gaps.
    //!\brief Stores the position (incl. gaps) where the last (consecutive) gap that is still before the current
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

        if (ungapped_view_pos != static_cast<int64_t>(host->ungapped_view.size()) && pos >= left_gap_end
            && (anchor_set_it == host->anchors.end() || pos < anchor_set_it->first))
            is_at_gap = false;
        else
            is_at_gap = true;
    }

public:
    /*!\name Associated types
     * \{
     */
    //!\brief The difference type.
    using difference_type = typename gap_decorator::difference_type;
    //!\brief The value type.
    using value_type = typename gap_decorator::value_type;
    //!\brief The reference type.
    using reference = typename gap_decorator::const_reference;
    //!\brief The pointer type.
    using pointer = value_type *;
    //!\brief The iterator category.
    using iterator_category = std::bidirectional_iterator_tag;
    //!\}

    /*!\name Constructors, destructor and assignment
     * \{
     */
    basic_iterator() = default;                                   //!< Defaulted.
    basic_iterator(basic_iterator const &) = default;             //!< Defaulted.
    basic_iterator & operator=(basic_iterator const &) = default; //!< Defaulted.
    basic_iterator(basic_iterator &&) = default;                  //!< Defaulted.
    basic_iterator & operator=(basic_iterator &&) = default;      //!< Defaulted.
    ~basic_iterator() = default;                                  //!< Defaulted.

    //!\brief Construct from seqan3::gap_decorator and initialising to first position.
    explicit basic_iterator(gap_decorator const & host_) : host(&host_), anchor_set_it{host_.anchors.begin()}
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
    basic_iterator(gap_decorator const & host_, typename gap_decorator::size_type const pos_) : host(&host_)
    {
        jump(pos_); // random access to pos
    }
    //!\}

    /*!\name Arithmetic operators
     * \{
     */
    //!\brief Increments iterator.
    basic_iterator & operator++()
    {
        assert(host); // host is set
        ++pos;

        if (pos < left_gap_end) // we stay within the preceding gap stretch
            return *this;

        if (anchor_set_it == host->anchors.end() || pos < anchor_set_it->first)
        { // proceed within the view since we are right of the previous gap but didn't arrive at the right gap yet
            ++ungapped_view_pos;
            if (ungapped_view_pos != static_cast<int64_t>(host->ungapped_view.size()))
                is_at_gap = false;
        }
        else
        { // we arrived at the right gap and have to update the variables. ungapped_view_pos remains unchanged.
            left_gap_end = anchor_set_it->first + anchor_set_it->second
                         - ((anchor_set_it != host->anchors.begin()) ? (std::prev(anchor_set_it))->second : 0);
            ++anchor_set_it;
            is_at_gap = true;

            if (left_gap_end == host->size()) // very last gap
                ++ungapped_view_pos;
        }

        return *this;
    }

    //!\brief Returns an incremented iterator copy.
    basic_iterator operator++(int)
    {
        basic_iterator cpy{*this};
        ++(*this);
        return cpy;
    }

    //!\brief Advances iterator by `skip` many positions.
    basic_iterator & operator+=(difference_type const skip)
    {
        this->jump(this->pos + skip);
        return *this;
    }

    //!\brief Returns an iterator copy advanced by `skip` many positions.
    basic_iterator operator+(difference_type const skip) const
    {
        return basic_iterator{*(this->host), this->pos + skip};
    }

    //!\brief Returns an iterator copy advanced by `skip` many positions.
    friend basic_iterator operator+(difference_type const skip, basic_iterator const & it)
    {
        return it + skip;
    }

    //!\brief Decrements iterator.
    basic_iterator & operator--()
    {
        assert(host); // host is set
        --pos;

        if (pos < left_gap_end)
        { // there was no gap before but we arrive at the left gap and have to update the variables.
            (anchor_set_it != host->anchors.begin()) ? --anchor_set_it : anchor_set_it;

            if (anchor_set_it != host->anchors.begin())
            {
                auto prev = std::prev(anchor_set_it);
                left_gap_end =
                    prev->first + prev->second - ((prev != host->anchors.begin()) ? std::prev(prev)->second : 0);
            }
            else // [[unlikely]]
            {
                left_gap_end = 0;
            }
            is_at_gap = true;
        }
        else if (anchor_set_it == host->anchors.end() || pos < anchor_set_it->first)
        { // we are neither at the left nor right gap
            --ungapped_view_pos;
            is_at_gap = false;
        }
        // else -> no op (we are still within the right gap stretch)

        return *this;
    }

    //!\brief Returns a decremented iterator copy.
    basic_iterator operator--(int)
    {
        basic_iterator cpy{*this};
        --(*this);
        return cpy;
    }

    //!\brief Advances iterator by `skip` many positions.
    basic_iterator & operator-=(difference_type const skip)
    {
        this->jump(this->pos - skip);
        return *this;
    }

    //!\brief Returns an iterator copy advanced by `skip` many positions.
    basic_iterator operator-(difference_type const skip) const
    {
        return basic_iterator{*(this->host), this->pos - skip};
    }

    //!\brief Returns an iterator copy advanced by `skip` many positions.
    friend basic_iterator operator-(difference_type const skip, basic_iterator const & it)
    {
        return it - skip;
    }

    //!\brief Returns the distance between two iterators.
    difference_type operator-(basic_iterator const lhs) const noexcept
    {
        return static_cast<difference_type>(this->pos - lhs.pos);
    }
    //!\}

    /*!\name Element access
     * \{
     */
    //!\brief Dereference operator returns a copy of the element currently pointed at.
    reference operator*() const
    {
        return (is_at_gap) ? reference{gap{}} : reference{host->ungapped_view[ungapped_view_pos]};
    }

    //!\brief Return underlying container value currently pointed at.
    reference operator[](difference_type const n) const
    {
        return *(*this + n);
    }
    //!\}

    /*!\name Comparison operators
     * \brief Compares iterators by virtual position.
     * \{
     */

    //!\brief Checks whether `*this` is equal to `rhs`.
    friend bool operator==(basic_iterator const & lhs, basic_iterator const & rhs) noexcept
    {
        return lhs.pos == rhs.pos;
    }

    //!\brief Checks whether `*this` is not equal to `rhs`.
    friend bool operator!=(basic_iterator const & lhs, basic_iterator const & rhs) noexcept
    {
        return lhs.pos != rhs.pos;
    }

    //!\brief Checks whether `*this` is less than `rhs`.
    friend bool operator<(basic_iterator const & lhs, basic_iterator const & rhs) noexcept
    {
        return lhs.pos < rhs.pos;
    }

    //!\brief Checks whether `*this` is greater than `rhs`.
    friend bool operator>(basic_iterator const & lhs, basic_iterator const & rhs) noexcept
    {
        return lhs.pos > rhs.pos;
    }

    //!\brief Checks whether `*this` is less than or equal to `rhs`.
    friend bool operator<=(basic_iterator const & lhs, basic_iterator const & rhs) noexcept
    {
        return lhs.pos <= rhs.pos;
    }

    //!\brief Checks whether `*this` is greater than or equal to `rhs`.
    friend bool operator>=(basic_iterator const & lhs, basic_iterator const & rhs) noexcept
    {
        return lhs.pos >= rhs.pos;
    }
    //!\}
};

} // namespace seqan3
