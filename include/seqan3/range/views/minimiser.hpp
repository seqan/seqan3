// -----------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2020, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2020, Knut Reinert & MPI für molekulare Genetik
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at: https://github.com/seqan/seqan3/blob/master/LICENSE.md
// -----------------------------------------------------------------------------------------------------

/*!\file
 * \author Mitra Darvish <mitra.darvish AT fu-berlin.de>
 * \brief Provides seqan3::views::minimiser.
 */

#pragma once

#include <deque>

#include <seqan3/range/concept.hpp>
#include <seqan3/range/views/detail.hpp>

namespace seqan3::detail
{
// ---------------------------------------------------------------------------------------------------------------------
// minimiser_view class
// ---------------------------------------------------------------------------------------------------------------------

/*!\brief The type returned by seqan3::views::minimiser.
 * \tparam urng_t The type of the underlying ranges, must model std::forward_range, the reference type must model
 *                std::totally_ordered. The typical use case is that the reference type is the result of
 *                seqan3::kmer_hash.
 * \implements std::ranges::view
 * \implements std::ranges::random_access_range
 * \implements std::ranges::sized_range
 * \ingroup views
 *
 * \details
 *
 * Note that most members of this class are generated by ranges::view_interface which is not yet documented here.
 */
template <std::ranges::view urng_t>
class minimiser_view : public std::ranges::view_interface<minimiser_view<urng_t>>
{
private:
    static_assert(std::ranges::forward_range<urng_t const>, "The minimiser_view only works on forward_ranges.");
    static_assert(std::totally_ordered<reference_t<urng_t>>, "The reference type of the underlying range must model "
                                                             "std::totally_ordered.");
    //!\brief The underlying range.
    urng_t urange;

    //!\brief The number of elements in one window.
    uint32_t num_w_elems;

    template <typename rng_t>
    class window_iterator;

public:
    /*!\name Constructors, destructor and assignment
     * \{
     */
    minimiser_view()                                       = default; //!< Defaulted.
    minimiser_view(minimiser_view const & rhs)             = default; //!< Defaulted.
    minimiser_view(minimiser_view && rhs)                  = default; //!< Defaulted.
    minimiser_view & operator=(minimiser_view const & rhs) = default; //!< Defaulted.
    minimiser_view & operator=(minimiser_view && rhs)      = default; //!< Defaulted.
    ~minimiser_view()                                      = default; //!< Defaulted.

    //!\brief Construct from a view and a given number of elements in one window.
    minimiser_view(urng_t urange_, uint32_t const & w_) : urange{std::move(urange_)}, num_w_elems{w_}{}

    //!\brief Construct from a non-view that can be view-wrapped and a given number of elements in one window.
    template <typename rng_t>
    //!\cond
     requires std::ranges::viewable_range<rng_t> &&
              std::constructible_from<urng_t, ranges::ref_view<std::remove_reference_t<rng_t>>>
    //!\endcond
    minimiser_view(rng_t && urange_, uint32_t const & w_) :
        urange{std::views::all(std::forward<rng_t>(urange_))}, num_w_elems{w_} {}
    //!\}

    /*!\name Iterators
     * \{
     */
    /*!\brief Returns an iterator to the first element of the range.
     * \returns Iterator to the first element.
     *
     * \details
     *
     * ### Complexity
     *
     * Constant.
     *
     * ### Exceptions
     *
     * No-throw guarantee.
     */
    auto begin() noexcept
    {
        return window_iterator<urng_t>{std::ranges::begin(urange), std::ranges::end(urange), num_w_elems};
    }

    //!\copydoc begin()
    auto begin() const noexcept
    //!\cond
        requires seqan3::const_iterable_range<urng_t>
    //!\endcond
    {
        return window_iterator<urng_t const>{std::ranges::begin(urange), std::ranges::end(urange), num_w_elems};
    }

    //!\copydoc begin()
    auto cbegin() const noexcept
    //!\cond
        requires seqan3::const_iterable_range<urng_t>
    //!\endcond
    {
        return begin();
    }

    /*!\brief Returns an iterator to the element following the last element of the range.
     * \returns Iterator to the end.
     *
     * \details
     *
     * This element acts as a placeholder; attempting to dereference it results in undefined behaviour.
     *
     * ### Complexity
     *
     * Constant.
     *
     * ### Exceptions
     *
     * No-throw guarantee.
     */
    auto end() noexcept
    {
        return std::ranges::end(urange);
    }

    //!\copydoc end()
    auto end() const noexcept
    //!\cond
        requires seqan3::const_iterable_range<urng_t>
    //!\endcond
    {
        return std::ranges::end(urange);
    }

    //!\copydoc end()
    auto cend() const noexcept
    //!\cond
        requires seqan3::const_iterable_range<urng_t>
    //!\endcond
    {
        return end();
    }
    //!\}
};

//!\brief Iterator for calculating minimisers.
template <std::ranges::view urng_t>
template <typename rng_t>
class minimiser_view<urng_t>::window_iterator
{
private:
    //!\brief The iterator type of the underlying range.
    using it_t = std::ranges::iterator_t<rng_t>;
    //!\brief The sentinel type of the underlying range.
    using sentinel_t = std::ranges::sentinel_t<rng_t>;

public:
    /*!\name Associated types
     * \{
     */
    //!\brief Type for distances between iterators.
    using difference_type = typename std::iter_difference_t<it_t>;
    //!\brief Value type of this iterator.
    using value_type = size_t;
    //!\brief The pointer type.
    using pointer = void;
    //!\brief Reference to `value_type`.
    using reference = value_type;
    //!\brief Tag this class as forward iterator.
    using iterator_category = std::forward_iterator_tag;
    //!\brief Tag this class depending on which concept `it_t` models.
    using iterator_concept = std::conditional_t<std::contiguous_iterator<it_t>,
                                                typename std::forward_iterator_tag,
                                                seqan3::iterator_tag_t<it_t>>;
    //!\}

    /*!\name Constructors, destructor and assignment
     * \{
     */
    constexpr window_iterator()                                    = default; //!< Defaulted.
    constexpr window_iterator(window_iterator const &)             = default; //!< Defaulted.
    constexpr window_iterator(window_iterator &&)                  = default; //!< Defaulted.
    constexpr window_iterator & operator=(window_iterator const &) = default; //!< Defaulted.
    constexpr window_iterator & operator=(window_iterator &&)      = default; //!< Defaulted.
    ~window_iterator()                                             = default; //!< Defaulted.

    /*!\brief Construct from begin and end iterators of a given range over std::totally_ordered elements, and the number
    *                      of elements per window.
    * /param[in] it_start Iterator pointing to the first position of the std::totally_ordered range.
    * /param[in] it_end   Iterator pointing to the last position of the std::totally_ordered range.
    * /param[in] w        The number of elements in one window.
    *
    * \details
    *
    * Looks at the number of elements per window, returns the smallest as minimiser and shifts then by one to repeat
    * this action. If a minimiser in consecutive windows is the same, it is returned only once.
    *
    */
    window_iterator(it_t it_start, sentinel_t it_end, uint32_t w) :
                    urange_end{it_end}, window_right{it_start}, num_w_elems{w}
    {
        if (num_w_elems > std::ranges::distance(window_right, urange_end))
            num_w_elems = std::ranges::distance(window_right, urange_end);
        if (window_right != urange_end)
            get_minimiser();
    }
    //!\}

    //!\anchor window_iterator_comparison
    //!\name Comparison operators
    //!\{

    //!\brief Compare to the sentinel of the underlying range.
    friend bool operator==(window_iterator const & lhs, sentinel_t const & rhs) noexcept
    {
        return lhs.window_right == rhs;
    }

    //!\brief Compare to iterator on underlying range.
    friend bool operator==(sentinel_t const & lhs, window_iterator const & rhs) noexcept
    {
        return lhs == rhs.window_right;
    }

    //!\brief Compare to the sentinel of the underlying range.
    friend bool operator==(window_iterator const & lhs, window_iterator const & rhs) noexcept
    {
        return std::tie(lhs.window_right, lhs.num_w_elems) == std::tie(rhs.window_right, rhs.num_w_elems);
    }

    //!\brief Compare to the sentinel of the underlying range.
    friend bool operator!=(window_iterator const & lhs, sentinel_t const & rhs) noexcept
    {
        return !(lhs == rhs);
    }

    //!\brief Compare to the sentinel of the underlying range.
    friend bool operator!=(sentinel_t const & lhs, window_iterator const & rhs) noexcept
    {
        return !(lhs == rhs);
    }

    //!\brief Compare to another window_iterator.
    friend bool operator!=(window_iterator const & lhs, window_iterator const & rhs) noexcept
    {
        return !(lhs == rhs);
    }

    //!\brief Compare to another window_iterator.
    friend bool operator<(window_iterator const & lhs, window_iterator const & rhs) noexcept
    {
        return (lhs.window_right < rhs.window_right) && (lhs.num_w_elems < rhs.num_w_elems);
    }

    //!\brief Compare to another window_iterator.
    friend bool operator>(window_iterator const & lhs, window_iterator const & rhs) noexcept
    {
        return (lhs.window_right > rhs.window_right) && (lhs.num_w_elems > rhs.num_w_elems);
    }

    //!\brief Compare to another window_iterator.
    friend bool operator<=(window_iterator const & lhs, window_iterator const & rhs) noexcept
    {
        return (lhs.window_right <= rhs.window_right) && (lhs.num_w_elems <= rhs.num_w_elems);
    }

    //!\brief Compare to another window_iterator.
    friend bool operator>=(window_iterator const & lhs, window_iterator const & rhs) noexcept
    {
        return (lhs.window_right >= rhs.window_right) && (lhs.num_w_elems >= rhs.num_w_elems);
    }
    //!\}

    //!\brief Pre-increment.
    window_iterator & operator++() noexcept
    {
        get_minimiser();
        return *this;
    }

    //!\brief Post-increment.
    window_iterator operator++(int) noexcept
    {
        window_iterator tmp{*this};
        get_minimiser();
        return tmp;
    }

    //!\brief Return the minimiser.
    value_type operator*() const noexcept
    {
        return minimiser_value;
    }

private:
    //!brief Iterator to last element in range.
    sentinel_t urange_end;

    //!\brief The minimiser value.
    size_t minimiser_value{0};

    //!\brief Iterator to the rightmost element of one window.
    it_t window_right;

    //!\brief The number of elements in one window.
    uint32_t num_w_elems;

    //!\brief Stored values per window. It is necessary to store them, because a shift can remove the current minimiser.
    std::deque<uint64_t> window_values;

    //!\brief Increments iterator by 1.
    void get_minimiser()
    {
        if (window_values.size() == 0)
            window_first();
        else // Call next_minimiser until minimiser value changed or end of the underlying range is reached.
            while(!next_minimiser()){}
    }

    //!\brief Calculates minimisers for the first window.
    void window_first()
    {
        for (uint32_t i = 0; (i < num_w_elems - 1) ; i++)
        {
            window_values.push_back(*window_right);
            std::ranges::advance(window_right,  1);
        }
        window_values.push_back(*window_right);
        minimiser_value = *(std::min_element(std::begin(window_values), std::end(window_values)));
    }

    //!\brief Calculates the next minimiser value.
    // For the following windows, we remove the first window element (is now not in window_values) and add the new
    // element that results from the window shifting.
    bool next_minimiser()
    {
        std::ranges::advance(window_right, 1);
        if (window_right == urange_end)
            return true;

        uint64_t new_element = *window_right;
        if (minimiser_value == *(std::begin(window_values)))
        {
            window_values.pop_front();
            if (!window_values.empty())
            {
                window_values.push_back(new_element);
                minimiser_value = *(std::min_element(std::begin(window_values), std::end(window_values)));
                return true;
            }

        }
        else
        {
            window_values.pop_front();
        }

        window_values.push_back(new_element);

        if (new_element < minimiser_value)
        {
            minimiser_value = new_element;
            return true;
        }
        else if (window_values.size() == 1)
        {
            return true;
        }

        return false;
    }
};

//!\brief A deduction guide for the view class template.
template <std::ranges::viewable_range rng_t>
minimiser_view(rng_t &&, uint32_t const & num_w_elems) -> minimiser_view<std::ranges::all_view<rng_t>>;


// ---------------------------------------------------------------------------------------------------------------------
// minimiser_fn (adaptor definition)
// ---------------------------------------------------------------------------------------------------------------------

//![adaptor_def]
//!\brief views::minimiser's range adaptor object type (non-closure).
struct minimiser_fn
{
    //!\brief Store the number of elements in one window and return a range adaptor closure object.
    constexpr auto operator()(uint32_t const & num_w_elems) const
    {
        return adaptor_from_functor{*this, num_w_elems};
    }

    /*!\brief                   Call the view's constructor with two arguments: the underlying view and an integer
     *                          indicating how many elements one window contains.
     * \param[in] urange        The input range to process. Must model std::ranges::viewable_range and
     *                          std::ranges::forward_range.
     * \param[in] num_w_elems   The number of elements in one window.
     * \returns                 A range of converted elements.
     */
    template <std::ranges::range urng_t>
    constexpr auto operator()(urng_t && urange, uint32_t const & num_w_elems) const
    {
        static_assert(std::ranges::viewable_range<urng_t>,
            "The range parameter to views::minimiser cannot be a temporary of a non-view range.");
        static_assert(std::ranges::forward_range<urng_t>,
            "The range parameter to views::minimiser must model std::ranges::forward_range.");

        return minimiser_view{urange, num_w_elems};
    }
};
//![adaptor_def]

} // namespace seqan3::detail

namespace seqan3::views
{

/*!\name Alphabet related views
 * \{
 */

/*!\brief                   Computes minimisers for a range of integral numbers. A minimiser is the smallest value in a
 *                          window.
 * \tparam urng_t           The type of the range being processed. See below for requirements. [template parameter is
 *                          omitted in pipe notation]
 * \param[in] urange        The range being processed. [parameter is omitted in pipe notation]
 * \param[in] num_w_elems   The number of elements in one window.
 * \returns                 A range of std::totally_ordered where each value is the minimal value for one window.
 *                          See below for the properties of the returned range.
 * \ingroup views
 *
 * \details
 *
 * A minimiser is the smallest value in a window. For example for the following list of hash values
 * [28, 100, 9, 23, 4, 1, 72, 37, 8] and 4 elements per window, the minimiser values are [9,4,1].
 *
 *
 * ### View properties
 *
 * | Concepts and traits              | `urng_t` (underlying range type)   | `rrng_t` (returned range type)   |
 * |----------------------------------|:----------------------------------:|:--------------------------------:|
 * | std::ranges::input_range         | *required*                         | *preserved*                      |
 * | std::ranges::forward_range       | *required*                         | *preserved*                      |
 * | std::ranges::bidirectional_range |                                    | *lost*                           |
 * | std::ranges::random_access_range |                                    | *lost*                           |
 * | std::ranges::contiguous_range    |                                    | *lost*                           |
 * |                                  |                                    |                                  |
 * | std::ranges::viewable_range      | *required*                         | *guaranteed*                     |
 * | std::ranges::view                |                                    | *guaranteed*                     |
 * | std::ranges::sized_range         |                                    | *lost*                           |
 * | std::ranges::common_range        |                                    | *lost*                           |
 * | std::ranges::output_range        |                                    | *lost*                           |
 * | seqan3::const_iterable_range     |                                    | *preserved*                      |
 * |                                  |                                    |                                  |
 * | std::ranges::range_reference_t   | std::totally_ordered               | std::totally_ordered             |
 *
 * See the \link views views submodule documentation \endlink for detailed descriptions of the view properties.
 *
 * ### Example
 *
 * \include test/snippet/range/views/minimiser.cpp
 *
 * \hideinitializer
 */
inline constexpr auto minimiser = detail::minimiser_fn{};

//!\}

} // namespace seqan3::views
