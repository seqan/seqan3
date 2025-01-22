// SPDX-FileCopyrightText: 2006-2025 Knut Reinert & Freie Universität Berlin
// SPDX-FileCopyrightText: 2016-2025 Knut Reinert & MPI für molekulare Genetik
// SPDX-License-Identifier: BSD-3-Clause

/*!\file
 * \author Rene Rahn <rene.rahn AT fu-berlin.de>
 * \brief Provides seqan3::views::single_pass_input
 */

#pragma once

#include <cassert>
#include <concepts>
#include <iterator>
#include <memory>
#include <ranges>
#include <type_traits>

#include <seqan3/core/detail/iterator_traits.hpp>
#include <seqan3/core/range/detail/adaptor_for_view_without_args.hpp>

//-----------------------------------------------------------------------------
// Implementation of single pass input view.
//-----------------------------------------------------------------------------

namespace seqan3::detail
{

//!\brief Forward declaration.
template <typename view_t>
class basic_iterator;

/*!\brief Adds single_pass_input behavior to the underlying range.
 * \tparam urng_t The underlying range type.
 * \implements std::ranges::input_range
 * \ingroup utility_views
 */
//![view_def]
template <std::ranges::view urng_t>
class single_pass_input_view : public std::ranges::view_interface<single_pass_input_view<urng_t>>
{
    //![view_def]

private:
    //!\brief The iterator type for the underlying range.
    using urng_iterator_type = std::ranges::iterator_t<urng_t>;

    //!\brief Friend declaration for seqan3::detail::basic_iterator.
    template <typename view_t>
    friend class basic_iterator;

    //!\brief An internal state to capture the underlying range and a cached iterator.
    struct state
    {
        //!\brief The underlying range.
        urng_t urng;
        //!\brief The cached iterator of the underlying range.
        urng_iterator_type cached_urng_iter = std::ranges::begin(urng);
    };

    //!\brief Manages the internal state.
    std::shared_ptr<state> state_ptr{};

    /*!\name Member types
     * \{
     */
    //!\brief Iterator type.
    using iterator = basic_iterator<single_pass_input_view>;
    //!\brief The sentinel type.
    using sentinel = std::ranges::sentinel_t<urng_t>;
    //\}

public:
    /*!\name Constructor, destructor, and assignment.
     * \{
     * \brief All standard functions are explicitly set to default.
     */
    //!\brief Default default-constructor.
    constexpr single_pass_input_view() = default;
    //!\brief Default copy-constructor.
    constexpr single_pass_input_view(single_pass_input_view const &) = default;
    //!\brief Default move-constructor.
    constexpr single_pass_input_view(single_pass_input_view &&) = default;
    //!\brief Default copy-assignment.
    constexpr single_pass_input_view & operator=(single_pass_input_view const &) = default;
    //!\brief Default move-assignment
    constexpr single_pass_input_view & operator=(single_pass_input_view &&) = default;
    //!\brief Default destructor.
    ~single_pass_input_view() = default;

    //!\brief Construction from the underlying view.
    explicit constexpr single_pass_input_view(urng_t _urng) : state_ptr{new state{std::move(_urng)}}
    {}

    //!\brief Construction from std::ranges::viewable_range.
    template <typename other_urng_t>
        requires (!std::same_as<std::remove_cvref_t<other_urng_t>, single_pass_input_view>
                  && std::ranges::viewable_range<other_urng_t>
                  && // Must come after self type check to avoid conflicts with the move constructor.
                  std::constructible_from<urng_t, std::ranges::ref_view<std::remove_reference_t<other_urng_t>>>)
    explicit constexpr single_pass_input_view(other_urng_t && _urng) : single_pass_input_view{std::views::all(_urng)}
    {}
    //!\}

    /*!\name Iterators
     * \{
     */
    /*!\brief Returns an iterator to the current begin of the underlying range.
     *
     * \details
     *
     * Subsequent calls to begin will result in different positions if the iterator was incremented
     * between the calls.
     */
    iterator begin()
    {
        return {*this};
    }

    //!\brief Const version of begin is deleted, since the underlying view_state must be mutable.
    iterator begin() const = delete;

    //!\brief Returns a sentinel.
    sentinel end()
    {
        return {std::ranges::end(state_ptr->urng)};
    }

    //!\brief Const version of end is deleted, since the underlying view_state must be mutable.
    sentinel end() const = delete;
    //!\}
};

/*!\name Deduction guide.
 * \relates seqan3::detail::single_pass_input_view
 * \{
 */

//!\brief Deduces the single_pass_input_view from the underlying range if it is a std::ranges::viewable_range.
template <std::ranges::viewable_range urng_t>
single_pass_input_view(urng_t &&) -> single_pass_input_view<std::views::all_t<urng_t>>;
//!\}
} // namespace seqan3::detail

//-----------------------------------------------------------------------------
// Iterator for single pass input view.
//-----------------------------------------------------------------------------

namespace seqan3::detail
{
/*!\brief An input_iterator over the associated range.
 * \implements std::input_iterator
 * \ingroup utility_views
 * \tparam view_type The type of the associated type.
 *
 * This iterator reduces every iterator type of the associated view to an single pass input iterator.
 */
template <typename view_type>
class basic_iterator<single_pass_input_view<view_type>>
{
    //!\brief The pointer to the associated view.
    using base_iterator_type = typename single_pass_input_view<view_type>::urng_iterator_type;
    //!\brief The sentinel type to compare to.
    using sentinel_type = typename single_pass_input_view<view_type>::sentinel;

    //!\brief The pointer to the associated view.
    single_pass_input_view<view_type> * view_ptr{};

    //!\brief Friend declaration to give seqan3::detail::single_pass_input_sentinel access to members of this class.
    template <typename input_view_type>
    friend class basic_iterator;

    //!\brief Test that the sentinel fulfills the std::sentinel_for for the underlying iterator.
    static_assert(std::sentinel_for<sentinel_type, base_iterator_type>);

public:
    /*!\name Associated types
     * \{
     */
    //!\brief Difference type.
    using difference_type = std::iter_difference_t<base_iterator_type>;
    //!\brief Value type.
    using value_type = std::iter_value_t<base_iterator_type>;
    //!\brief Pointer type.
    using pointer = detail::iter_pointer_t<base_iterator_type>;
    //!\brief Reference type.
    using reference = std::iter_reference_t<base_iterator_type>;
    //!\brief Iterator category.
    using iterator_category = std::input_iterator_tag;
    //!\}

    /*!\name Construction, destruction and assignment
     * \{
     */
    //!\brief Default construction.
    basic_iterator() = default;
    //!\brief Copy construction.
    constexpr basic_iterator(basic_iterator const & rhs) = default;
    //!\brief Move construction.
    constexpr basic_iterator(basic_iterator && rhs) = default;
    //!\brief Copy assignment.
    constexpr basic_iterator & operator=(basic_iterator const & rhs) = default;
    //!\brief Move assignment.
    constexpr basic_iterator & operator=(basic_iterator && rhs) = default;
    //!\brief Destruction.
    ~basic_iterator() = default;

    //!\brief Constructing from the underlying seqan3::single_pass_input_view.
    basic_iterator(single_pass_input_view<view_type> & view) noexcept : view_ptr{&view}
    {}
    //!\}

    /*!\name Access operations
     * \{
     */
    //!\brief Dereferences the cached iterator.
    reference operator*() const noexcept
    {
        return *cached();
    }

    //!\brief Returns pointer to the pointed-to object.
    pointer operator->() const noexcept
        requires (!std::is_void_v<pointer>)
    {
        return std::addressof(*cached());
    }
    //!\}

    /*!\name Iterator operations
     * \{
     */
    //!\brief Pre-increment.
    basic_iterator & operator++() noexcept
    {
        ++cached();
        return *this;
    }

    //!\brief Post-increment.
    void operator++(int) noexcept
    {
        // this post-increment can't be an std::output_iterator, because it would require that `*it++ = value` must have
        // the same semantic as `*i = value; ++i`, but it actually has the following `++i; *i = value` semantic, due to
        // the centralised storage of the underlying_iterator in the view where each copy of a basic_iterator points to
        // the same centralised state.
        ++(*this);
    }
    //!\}

    /*!\name Comparison operators
     * \{
     */
    //!\brief Compares for equality with sentinel.
    constexpr bool operator==(sentinel_type const & s) const noexcept
    {
        return cached() == s;
    }

    //!\copydoc operator==
    friend constexpr bool operator==(sentinel_type const & s, basic_iterator const & rhs) noexcept
    {
        return rhs == s;
    }

    //!\brief Compares for inequality with sentinel.
    constexpr bool operator!=(sentinel_type const & rhs) const noexcept
    {
        return !(*this == rhs);
    }

    //!\copydoc operator!=
    friend constexpr bool operator!=(sentinel_type const & s, basic_iterator const & rhs) noexcept
    {
        return rhs != s;
    }
    //!\}

protected:
    //!\privatesection
    //!\brief Gives access to the cached iterator.
    base_iterator_type & cached() const noexcept
    {
        assert(view_ptr != nullptr);
        assert(view_ptr->state_ptr != nullptr);
        return view_ptr->state_ptr->cached_urng_iter;
    }
};
} // namespace seqan3::detail

//-----------------------------------------------------------------------------
// View shortcut for functor.
//-----------------------------------------------------------------------------

//![adaptor_def]
namespace seqan3::views
{
/*!\brief               A view adapter that decays most of the range properties and adds single pass behavior.
 * \tparam urng_t       The type of the range being processed. See below for requirements.
 * \param[in] urange    The range being processed.
 * \returns             A range with single pass input behavior. See below for the properties of the returned range.
 * \ingroup utility_views
 *
 * \details
 *
 * \header_file{seqan3/utility/views/single_pass_input.hpp}
 *
 * This view adds single-pass semantics to any input view. This means, that `begin` always returns the iterator
 * to the current location in the underlying range after `k` elements have been already consumed and not to the begin
 * of the underlying range, i.e. it mirrors the behavior of an input stream.
 * Note, the view updates an internal state after moving the associated iterator.
 * Thus, the `const begin` and `const end` are explicitly deleted.
 *
 * ### View properties
 *
 *
 * | Concepts and traits              | `urng_t` (underlying range type)      | `rrng_t` (returned range type)                     |
 * |----------------------------------|:-------------------------------------:|:--------------------------------------------------:|
 * | std::ranges::input_range         | *required*                            | *preserved*                                        |
 * | std::ranges::forward_range       |                                       | *lost*                                             |
 * | std::ranges::bidirectional_range |                                       | *lost*                                             |
 * | std::ranges::random_access_range |                                       | *lost*                                             |
 * | std::ranges::contiguous_range    |                                       | *lost*                                             |
 * |                                  |                                       |                                                    |
 * | std::ranges::viewable_range      | *required*                            | *guaranteed*                                       |
 * | std::ranges::view                |                                       | *guaranteed*                                       |
 * | std::ranges::sized_range         |                                       | *lost*                                             |
 * | std::ranges::common_range        |                                       | *lost*                                             |
 * | std::ranges::output_range        |                                       | *lost*                                             |
 * | seqan3::const_iterable_range     |                                       | *lost*                                             |
 * |                                  |                                       |                                                    |
 * | std::ranges::range_reference_t   |                                       | std::ranges::range_reference_t<urng_t>             |
 *
 * See the \link views views submodule documentation \endlink for detailed descriptions of the view properties.
 *
 * ### Thread safety
 *
 * Concurrent access to this view, e.g. while iterating over it, is not thread-safe and must be protected externally.
 *
 * ### Example
 *
 * \include test/snippet/utility/views/single_pass_input.cpp
 * \hideinitializer
 */
inline constexpr auto single_pass_input = detail::adaptor_for_view_without_args<detail::single_pass_input_view>{};

} // namespace seqan3::views
//![adaptor_def]
