// -----------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2020, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2020, Knut Reinert & MPI für molekulare Genetik
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at: https://github.com/seqan/seqan3/blob/master/LICENSE.md
// -----------------------------------------------------------------------------------------------------

/*!\file
 * \brief Provides seqan3::detail::counted_simd_iterator and seqan3::views::iota_simd.
 * \author Rene Rahn <rene.rahn AT fu-berlin.de>
 */

#pragma once

#include <seqan3/std/ranges>

#include <seqan3/core/concept/core_language.hpp>
#include <seqan3/core/simd/concept.hpp>
#include <seqan3/core/simd/simd_algorithm.hpp>
#include <seqan3/core/simd/simd_traits.hpp>

namespace seqan3::detail
{

/*!\brief Implements a special version of a counted iterator over a simd vector.
 * \ingroup simd
 * \implements std::forward_iterator
 *
 * \tparam index_t The type of the index; must model seqan3::simd::simd_concept.
 *
 * \details
 *
 * Uses a simd count vector to increment the counted iterator. This seems to be in general faster than calling
 * seqan3::simd::fill when dereferencing the iterator, although this is just a constant and fast operation.
 */
template <simd_concept index_t>
class counted_simd_iterator
{
private:
    //!\brief The currently represented count.
    index_t count_simd{};
    //!\brief The count in scalar representation.
    size_t count_scalar{};

public:
    /*!\name Associated types
     * \{
     */
    //!\brief The value type.
    using value_type = index_t;
    //!\brief The reference type.
    using reference = value_type;
    //!\brief The pointer type.
    using pointer = void;
    //!\brief The difference type.
    using difference_type = std::ptrdiff_t;
    //!\brief The iterator category.
    using iterator_category = std::forward_iterator_tag;
    //!\}

    /*!\name Constructor, assignment and destructor
     * \{
     */
    counted_simd_iterator() noexcept = default; //!< Defaulted.
    counted_simd_iterator(counted_simd_iterator const &) noexcept = default; //!< Defaulted.
    counted_simd_iterator(counted_simd_iterator &&) noexcept = default; //!< Defaulted.
    counted_simd_iterator & operator=(counted_simd_iterator const &) noexcept = default; //!< Defaulted.
    counted_simd_iterator & operator=(counted_simd_iterator &&) noexcept = default; //!< Defaulted.
    ~counted_simd_iterator() = default; //!< Defaulted.

    /*!\brief Constructs and initialises the iterator with the given index.
     *
     * \tparam scalar_index_t The scalar index type; must model seqan3::arithmetic.
     *
     * \param[in] scalar_index The scalar index the iterator shall represent.
     */
    template <arithmetic scalar_index_t>
    explicit counted_simd_iterator(scalar_index_t const scalar_index) noexcept :
        count_simd{simd::fill<index_t>(scalar_index)},
        count_scalar{static_cast<size_t>(scalar_index)}
    {}
    //!\}

    /*!\name Element access
     * \{
     */

    //!\brief Access the pointed-to matrix coordinate column.
    reference operator*() const
    {
        return count_simd;
    }
    //!\}

    /*!\name Arithmetic operators
     * \{
     */

    //!\brief Increments the iterator.
    counted_simd_iterator & operator++()
    {
        count_simd += seqan3::simd::fill<index_t>(1);
        ++count_scalar;
        return *this;
    }

    //!\brief Increments the iterator and returns the iterator pointing to the previous index.
    counted_simd_iterator operator++(int)
    {
        counted_simd_iterator tmp{*this};
        ++(*this);
        return tmp;
    }

    //!\brief Returns the distance between two iterators.
    difference_type operator-(counted_simd_iterator const & rhs) const
    {
        return count_scalar - rhs.count_scalar;
    }
    //!\}

    /*!\name Comparison operators
     * \{
     */

    //!\brief Tests whether `lhs == rhs`.
    friend bool operator==(counted_simd_iterator const & lhs, counted_simd_iterator const & rhs) noexcept
    {
        return lhs.count_scalar == rhs.count_scalar;
    }

    //!\brief Tests whether `lhs != rhs`.
    friend bool operator!=(counted_simd_iterator const & lhs, counted_simd_iterator const & rhs) noexcept
    {
        return !(lhs == rhs);
    }
    //!\}
};

} // namespace seqan3::detail

namespace seqan3::views
{

/*!\brief A closed iota view over a simd vector.
 * \ingroup simd
 * \implements std::ranges::forward_range
 * \implements std::ranges::sized_range
 * \implements std::ranges::common_range
 * \implements std::ranges::borrowed_range
 *
 * \tparam index_t The represented index type; must model seqan3::simd::simd_concept.
 *
 * \details
 *
 * This view is a lightweight wrapper around a seqan3::detail::counted_simd_iterator pair.
 * The regular std::views::iota view cannot be used in combination with a simd vector type because the comparison of
 * two simd vectors does not return a bool but another simd vector type.
 *
 * ### View properties
 *
 * | Concepts and traits              | `rrng_t` (returned range type) |
 * |----------------------------------|:------------------------------:|
 * | std::ranges::input_range         | *guaranteed*                   |
 * | std::ranges::forward_range       | *guaranteed*                   |
 * | std::ranges::bidirectional_range | *lost*                         |
 * | std::ranges::random_access_range | *lost*                         |
 * | std::ranges::contiguous_range    | *lost*                         |
 * |                                  |                                |
 * | std::ranges::viewable_range      | *guaranteed*                   |
 * | std::ranges::view                | *guaranteed*                   |
 * | std::ranges::sized_range         | *guaranteed*                   |
 * | std::ranges::common_range        | *guaranteed*                   |
 * | std::ranges::output_range        | *lost*                         |
 * | std::ranges::borrowed_range      | *guaranteed*                   |
 * | seqan3::const_iterable_range     | *guaranteed*                   |
 * |                                  |                                |
 * | std::ranges::range_reference_t   | simd_t                         |
 *
 * For more details, see \ref views.
 *
 * ### Example
 *
 * \include test/snippet/core/simd/view_iota_simd.cpp
 */
template <simd_concept index_t>
class iota_simd : public std::ranges::view_interface<iota_simd<index_t>>
{
private:
    //!\brief The underlying scalar type.
    using scalar_type = typename simd_traits<index_t>::scalar_type;
    //!\brief The counted simd iterator type.
    using iterator_type = seqan3::detail::counted_simd_iterator<index_t>;

    //!\brief The begin index.
    scalar_type begin_index{};
    //!\brief The end index.
    scalar_type end_index{};

public:
    /*!\name Constructor, assignment and destructor
     * \{
     */
    iota_simd() noexcept = default; //!< Defaulted.
    iota_simd(iota_simd const &) noexcept = default; //!< Defaulted.
    iota_simd(iota_simd &&) noexcept = default; //!< Defaulted.
    iota_simd & operator=(iota_simd const &) noexcept = default; //!< Defaulted.
    iota_simd & operator=(iota_simd &&) noexcept = default; //!< Defaulted.
    ~iota_simd() = default; //!< Defaulted.

    /*!\brief Constructs the iota view from the given index pair.
     *
     * param[in] begin_index The first index of the iota view.
     * param[in] end_index The index behind the last one of the iota view.
     */
    iota_simd(scalar_type const begin_index, scalar_type const end_index) :
        begin_index{begin_index},
        end_index{end_index}
    {}
    //!\}

    /*!\name Iterators
     * \{
     */
    //!\brief Returns seqan3::detail::counted_simd_iterator pointing to the begin of the range.
    iterator_type begin() const noexcept
    {
        return iterator_type{begin_index};
    }

    //!\brief Returns seqan3::detail::counted_simd_iterator pointing to the end of the range.
    iterator_type end() const noexcept
    {
        return iterator_type{end_index};
    }
    //!\}
};

} // namespace seqan3::views

//!\cond
// Enable borrowed range.
// namespace std::ranges
// {
// template <typename index_t>
// inline constexpr bool enable_borrowed_range<seqan3::views::iota_simd<index_t>> = true;
// }
//!\endcond
