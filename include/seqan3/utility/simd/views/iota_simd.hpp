// SPDX-FileCopyrightText: 2006-2025 Knut Reinert & Freie Universität Berlin
// SPDX-FileCopyrightText: 2016-2025 Knut Reinert & MPI für molekulare Genetik
// SPDX-License-Identifier: BSD-3-Clause

/*!\file
 * \brief Provides seqan3::detail::counted_simd_iterator and seqan3::views::iota_simd.
 * \author Rene Rahn <rene.rahn AT fu-berlin.de>
 */

#pragma once

#include <ranges>

#include <seqan3/utility/concept.hpp>
#include <seqan3/utility/simd/algorithm.hpp>
#include <seqan3/utility/simd/concept.hpp>
#include <seqan3/utility/simd/simd_traits.hpp>

namespace seqan3::detail
{

/*!\brief Implements a special version of a counted iterator over a simd vector.
 * \ingroup utility_simd_views
 * \implements std::forward_iterator
 *
 * \tparam index_simd_t The type of the index; must model seqan3::simd::simd_concept.
 *
 * \details
 *
 * Uses a simd count vector to increment the counted iterator. This seems to be in general faster than calling
 * seqan3::simd::fill when dereferencing the iterator, although this is just a constant and fast operation.
 */
template <simd_concept index_simd_t>
class counted_simd_iterator
{
private:
    //!\brief The currently represented count.
    index_simd_t count_simd{};

public:
    /*!\name Associated types
     * \{
     */
    //!\brief The value type.
    using value_type = index_simd_t;
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
    counted_simd_iterator() = default;                                          //!< Defaulted.
    counted_simd_iterator(counted_simd_iterator const &) = default;             //!< Defaulted.
    counted_simd_iterator(counted_simd_iterator &&) = default;                  //!< Defaulted.
    counted_simd_iterator & operator=(counted_simd_iterator const &) = default; //!< Defaulted.
    counted_simd_iterator & operator=(counted_simd_iterator &&) = default;      //!< Defaulted.
    ~counted_simd_iterator() = default;                                         //!< Defaulted.

    /*!\brief Constructs and initialises the iterator with the given index.
     *
     * \tparam index_scalar_t The scalar index type; must model seqan3::arithmetic.
     *
     * \param[in] scalar_index The scalar index the iterator shall represent.
     */
    template <arithmetic index_scalar_t>
    explicit counted_simd_iterator(index_scalar_t const scalar_index) noexcept :
        count_simd{simd::fill<index_simd_t>(scalar_index)}
    {}
    //!\}

    /*!\name Element access
     * \{
     */

    //!\brief Return the current simd index.
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
        count_simd += seqan3::simd::fill<index_simd_t>(1);
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
        return count_simd[0] - rhs.count_simd[0];
    }
    //!\}

    /*!\name Comparison operators
     * \{
     */

    //!\brief Tests whether `lhs == rhs`.
    friend bool operator==(counted_simd_iterator const & lhs, counted_simd_iterator const & rhs) noexcept
    {
        return lhs.count_simd[0] == rhs.count_simd[0];
    }

    //!\brief Tests whether `lhs != rhs`.
    friend bool operator!=(counted_simd_iterator const & lhs, counted_simd_iterator const & rhs) noexcept
    {
        return !(lhs == rhs);
    }
    //!\}
};

/*!\brief The simd iota view.
 * \ingroup utility_simd_views
 * \tparam index_simd_t The represented index type; must model seqan3::simd::simd_concept.
 */
template <simd_concept index_simd_t>
class iota_simd_view : public std::ranges::view_interface<iota_simd_view<index_simd_t>>
{
private:
    //!\brief The underlying scalar type.
    using index_scalar_type = typename simd_traits<index_simd_t>::scalar_type;
    //!\brief The counted simd iterator type.
    using iterator_type = seqan3::detail::counted_simd_iterator<index_simd_t>;

    //!\brief The begin index.
    index_scalar_type begin_index{};
    //!\brief The end index.
    index_scalar_type end_index{};

public:
    /*!\name Constructor, assignment and destructor
     * \{
     */
    iota_simd_view() = default;                                   //!< Defaulted.
    iota_simd_view(iota_simd_view const &) noexcept = default;    //!< Defaulted.
    iota_simd_view(iota_simd_view &&) = default;                  //!< Defaulted.
    iota_simd_view & operator=(iota_simd_view const &) = default; //!< Defaulted.
    iota_simd_view & operator=(iota_simd_view &&) = default;      //!< Defaulted.
    ~iota_simd_view() = default;                                  //!< Defaulted.

    /*!\brief Constructs the iota view from the given index pair.
     *
     * param[in] begin_index The first index of the iota view.
     * param[in] end_index The index behind the last one of the iota view.
     */
    iota_simd_view(index_scalar_type const begin_index, index_scalar_type const end_index) :
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

/*!\brief The view adaptor returning the seqan3::detail::iota_simd_view.
 * \ingroup utility_simd_views
 * \tparam index_simd_t The represented index type; must model seqan3::simd::simd_concept.
 */
template <simd_concept index_simd_t>
struct iota_simd_view_fn
{
    //!\brief The underlying scalar type.
    using index_scalar_type = typename simd_traits<index_simd_t>::scalar_type;

    /*!\brief Returns a simd iota view over the given range.
     *
     * \param begin_index The begin index.
     * \param end_index The end index.
     *
     * \returns An instance of seqan3::detail::iota_simd_view.
     */
    constexpr auto operator()(index_scalar_type const begin_index, index_scalar_type const end_index) const
    {
        return iota_simd_view<index_simd_t>{begin_index, end_index};
    }
};
} // namespace seqan3::detail

namespace seqan3::views
{

/*!\brief An iota view over a simd vector.
 * \ingroup utility_simd_views
 * \implements std::ranges::forward_range
 * \implements std::ranges::sized_range
 * \implements std::ranges::common_range
 * \implements std::ranges::borrowed_range
 *
 * \tparam index_simd_t The represented index type; must model seqan3::simd::simd_concept.
 *
 * \details
 *
 * This view is an equivalent implementation to:
 *
 * \include test/snippet/utility/simd/views/iota_simd_transform.cpp
 *
 * However, benchmarks showed that increasing a simd vector is faster than constructing it every time
 * (up-to 2x speed-up). This speed-up justifies an own class that does this task more efficiently.
 *
 * This view is a lightweight wrapper around a seqan3::detail::counted_simd_iterator pair.
 * Note the regular std::views::iota view cannot be used with two simd types:
 *
 * ```cpp
 * std::views::iota(simd_begin, simd_end);
 * ```
 * because the return type of the comparison of two simd vector types is not convertible to `bool`.
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
 * | std::ranges::range_reference_t   | index_simd_t                   |
 *
 * This is a source view. For more details, see \ref views.
 *
 * ### Example
 *
 * \include test/snippet/utility/simd/views/iota_simd.cpp
 * \hideinitializer
 */
template <simd_concept index_simd_t>
inline constexpr detail::iota_simd_view_fn<index_simd_t> iota_simd{};

} // namespace seqan3::views

namespace std::ranges
{
//!\cond
template <seqan3::simd_concept index_simd_t>
inline constexpr bool enable_borrowed_range<seqan3::detail::iota_simd_view<index_simd_t>> = true;
//!\endcond
} // namespace std::ranges
