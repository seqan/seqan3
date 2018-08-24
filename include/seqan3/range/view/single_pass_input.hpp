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
 * \author Rene Rahn <rene.rahn AT fu-berlin.de>
 * \brief Provides seqan3::single_pass_input_view
 */

#pragma once

#include <seqan3/core/metafunction/all.hpp>
#include <seqan3/range/view/detail.hpp>
#include <seqan3/std/concepts>
#include <seqan3/std/iterator>
#include <seqan3/std/ranges>

//-----------------------------------------------------------------------------
// Implementation of single pass input view.
//-----------------------------------------------------------------------------

namespace seqan3::detail
{

//!\brief Forward declaration.
template <typename view_t>
class single_pass_input_iterator;

/*!\brief Adds single_pass_input behavior to the underlying range.
 * \tparam urng_t The underlying range type.
 * \implements std::ranges::InputRange
 * \ingroup view
 */
template <std::ranges::InputRange urng_t>
class single_pass_input_view : public detail::view_base
{
private:

    //!\brief The pure range type without any reference type.
    using pure_range_type    = std::remove_reference_t<urng_t>;
    //!\brief The iterator type for the underlying range.
    using urng_iterator_type = iterator_t<pure_range_type>;

    //!\brief Friend declaration for seqan3::detail::single_pass_input_iterator.
    template <typename view_t>
    friend class single_pass_input_iterator;

    //!\brief The internal view_state.
    struct view_state
    {
        //!\brief The underlying range.
        urng_t             urng;
        //!\brief The cached iterator of the underlying range.
        urng_iterator_type cached_urng_iter{};
    };

    //!\brief Shared pointer of the data model.
    std::shared_ptr<view_state> view_state_ptr{};

public:
    /*!\name Member types
     * \{
     */
    //!\brief Iterator type.
    using iterator          = single_pass_input_iterator<single_pass_input_view>;
    //!\brief Const iterator type is `void`, as iterating over this view as `const` is explicitly forbidden.
    using const_iterator    = void;
    //!\brief The sentinel type.
    using sentinel          = sentinel_t<pure_range_type>;
    //!\brief Value type.
    using value_type        = typename iterator::value_type;
    //!\brief Always returns immutable reference type, since single_pass_input cannot change the underlying values.
    using reference         = typename iterator::reference;
    //!\brief The const_reference type is `void`, as iterating over this view as `const` is explicitly forbidden.
    using const_reference   = void;
    //\}

    /*!\name Constructor, destructor, and assignment.
     * \{
     * \brief All standard functions are explicitly defaulted.
     */
    constexpr single_pass_input_view() = default;
    constexpr single_pass_input_view(single_pass_input_view const &) = default;
    constexpr single_pass_input_view(single_pass_input_view &&) = default;
    constexpr single_pass_input_view & operator=(single_pass_input_view const &) = default;
    constexpr single_pass_input_view & operator=(single_pass_input_view &&) = default;
    ~single_pass_input_view() = default;

    //!\brief Construction from the underlying view.
    single_pass_input_view(urng_t && urng) :
        view_state_ptr{new view_state{std::forward<urng_t>(urng), ranges::begin(urng)}}
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
        return iterator{*this};
    }

    //!\brief Const version of begin is deleted, since the underlying view_state must be mutable.
    const_iterator begin() const = delete;

    //!\copydoc single_pass_input_view::begin() const
    const_iterator cbegin() const = delete;

    //!\brief Returns a sentinel.
    sentinel end()
    {
        return {ranges::end(view_state_ptr->urng)};
    }

    //!\brief Const version of end is deleted, since the underlying view_state must be mutable.
    sentinel end() const = delete;

    //!\copydoc single_pass_input_view::end() const
    sentinel cend() const = delete;
    //!\}
};

/*!\name Deduction guide.
 * \relates seqan3::detail::single_pass_input_view
 * \{
 */
//!\brief Deduces the single_pass_input_view from the underlying range.
template <std::ranges::InputRange urng_t>
single_pass_input_view(urng_t &&) -> single_pass_input_view<urng_t>;
//!\}
} // seqan3::detail

//-----------------------------------------------------------------------------
// Iterator for single pass input view.
//-----------------------------------------------------------------------------

namespace seqan3::detail
{
/*!\brief An input_iterator over the associated range.
 * \implements std::InputIterator
 * \ingroup view
 * \tparam view_type The type of the associated type.
 *
 * This iterator reduces every iterator type of the associated view to an single pass input iterator.
 */
template <typename view_type>
class single_pass_input_iterator<single_pass_input_view<view_type>>
{
    //!\brief The pointer to the associated view.
    using base_iterator_type = typename single_pass_input_view<view_type>::urng_iterator_type;
    //!\brief The sentinel type to compare to.
    using sentinel_type      = typename single_pass_input_view<view_type>::sentinel;

    //!\brief The pointer to the associated view.
    single_pass_input_view<view_type> * view_ptr{};

    //!\ Friend declaration to give seqan3::detail::single_pass_input_sentinel access to members of this class.
    template <typename input_view_type>
    friend class single_pass_input_iterator;

    //!\brief Test that the sentinel fulfills the std::Sentinel for the underlying iterator.
    static_assert(std::Sentinel<sentinel_type, base_iterator_type>);

public:

    /*!\name Member types
     * \{
     */
    //!\brief Difference type.
    using difference_type   = difference_type_t<base_iterator_type>;
    //!\brief Value type.
    using value_type        = value_type_t<base_iterator_type>;
    //!\brief Pointer type.
    using pointer           = typename std::iterator_traits<base_iterator_type>::pointer;
    //!\brief Reference type.
    using reference         = reference_t<base_iterator_type>;
    //!\brief Iterator category.
    using iterator_category = std::input_iterator_tag;
    //!\}

    /*!\name Construction, destruction and assignment
     * \{
     */
    //!\brief Default construction.
    single_pass_input_iterator() = default;
    //!\brief Copy construction.
    constexpr single_pass_input_iterator(single_pass_input_iterator const & rhs) = default;
    //!\brief Move construction.
    constexpr single_pass_input_iterator(single_pass_input_iterator && rhs) = default;
    //!\brief Copy assignment.
    constexpr single_pass_input_iterator & operator=(single_pass_input_iterator const & rhs) = default;
    //!\brief Move assignment.
    constexpr single_pass_input_iterator & operator=(single_pass_input_iterator && rhs) = default;
    //!\brief Destruction.
    ~single_pass_input_iterator() = default;

    //!\brief Constructing from the underlying seqan3::single_pass_input_view.
    single_pass_input_iterator(single_pass_input_view<view_type> & view) noexcept : view_ptr{&view}
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
    //!\}

    /*!\name Iterator operations
     * \{
     */
    //!\brief Pre-increment.
    single_pass_input_iterator & operator++() noexcept
    {
        ++cached();
        return *this;
    }

    //!\brief Post-increment.
    void operator++(int) noexcept
    {
        ++(*this);
    }
    //!\}

    /*!\name Comparison operators
     * \{
     */
    //!\brief Compares iterator with sentinel.
    constexpr bool operator==(sentinel_type const & s) const noexcept
    {
        return cached() == s;
    }

    //!\copydoc operator==
    constexpr bool operator!=(sentinel_type const & rhs) const noexcept
    {
        return !(*this == rhs);
    }

    //!\copydoc operator==
    friend constexpr bool
    operator==(sentinel_type const & s,
               single_pass_input_iterator<single_pass_input_view<view_type>> const & rhs) noexcept
    {
        return rhs == s;
    }

    //!\copydoc operator==
    friend constexpr bool
    operator!=(sentinel_type const & s,
               single_pass_input_iterator<single_pass_input_view<view_type>> const & rhs) noexcept
    {
        return rhs != s;
    }
    //!\}

protected:
//!\privatesection
    //!\brief Gives access to the cached iterator.
    base_iterator_type & cached() const noexcept
    {
        return view_ptr->view_state_ptr->cached_urng_iter;
    }
};
}  // seqan3::detail

//-----------------------------------------------------------------------------
// View shortcut for functor.
//-----------------------------------------------------------------------------

namespace seqan3::view
{
/*!\name General purpose views
 * \{
 */

/*!\brief               A view adapter that decays most of the range properties and adds single pass behavior.
 * \tparam urng_t       The type of the range being processed. See below for requirements.
 * \param[in] urange    The range being processed.
 * \returns             A range with single pass input behavior. See below for the properties of the returned range.
 * \ingroup view
 *
 * \details
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
 * | range concepts and reference_t  | `urng_t` (underlying range type)      | `rrng_t` (returned range type)                     |
 * |---------------------------------|:-------------------------------------:|:--------------------------------------------------:|
 * | std::ranges::InputRange         | *required*                            | *preserved*                                        |
 * | std::ranges::ForwardRange       |                                       | *lost*                                             |
 * | std::ranges::BidirectionalRange |                                       | *lost*                                             |
 * | std::ranges::RandomAccessRange  |                                       | *lost*                                             |
 * | std::ranges::ContiguousRange    |                                       | *lost*                                             |
 * |                                 |                                       |                                                    |
 * | std::ranges::ViewableRange      | *required*                            | *guaranteed*                                       |
 * | std::ranges::View               |                                       | *guaranteed*                                       |
 * | std::ranges::SizedRange         |                                       | *lost*                                             |
 * | std::ranges::CommonRange        |                                       | *lost*                                             |
 * | std::ranges::OutputRange        |                                       | *preserved*                                        |
 * | seqan3::const_iterable_concept  |                                       | *lost*                                             |
 * |                                 |                                       |                                                    |
 * | seqan3::reference_t             |                                       | seqan3::reference_t<urng_t>                        |
 *
 * See the \link view view submodule documentation \endlink for detailed descriptions of the view properties.
 *
 * ### Thread safety
 *
 * Concurrent access to this view, e.g. while iterating over it, is not thread-safe and must be protected externally.
 *
 * ### Example
 *
 * ```cpp
 * std::string str{"hello"};
 * auto v = str | view::single_pass_input;
 * std::cout << *++v.begin() << std::endl;  // prints 'e'
 * std::cout << *++v.begin() << std::endl;  // prints 'l'
 * std::cout << *++v.begin() << std::endl;  // prints 'l'
 * std::cout << *++v.begin() << std::endl;  // prints 'o'
 * ```
 * \hideinitializer
 */
inline constexpr auto single_pass_input = detail::generic_pipable_view_adaptor<detail::single_pass_input_view>{};

//!\}
} // namespace seqan3::view
