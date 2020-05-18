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

    //!\brief The number of values in one window.
    uint32_t window_values_size;

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

    /*!\brief Construct from a view and a given number of values in one window.
    * \param[in] urange_   The input range to process. Must model std::ranges::viewable_range and
    *                      std::ranges::forward_range.
    * \param[in] w_        The number of values in one window.
    */
    minimiser_view(urng_t urange_, uint32_t const w_) : urange{std::move(urange_)}, window_values_size{w_}{}

    /*!\brief Construct from a non-view that can be view-wrapped and a given number of values in one window.
    * \param[in] urange_   The input range to process. Must model std::ranges::viewable_range and
    *                      std::ranges::forward_range.
    * \param[in] w_        The number of values in one window.
    */
    template <typename rng_t>
    //!\cond
     requires std::ranges::viewable_range<rng_t> &&
              std::constructible_from<urng_t, ranges::ref_view<std::remove_reference_t<rng_t>>>
    //!\endcond
    minimiser_view(rng_t && urange_, uint32_t const w_) :
    urange{std::views::all(std::forward<rng_t>(urange_))}, window_values_size{w_} {}
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
     * Strong exception guarantee.
     */
    auto begin()
    {
        return window_iterator<urng_t>{std::ranges::begin(urange), std::ranges::end(urange), window_values_size};
    }

    //!\copydoc begin()
    auto begin() const
    //!\cond
        requires seqan3::const_iterable_range<urng_t>
    //!\endcond
    {
        return window_iterator<urng_t const>{std::ranges::begin(urange), std::ranges::end(urange), window_values_size};
    }

    //!\copydoc begin()
    auto cbegin() const
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
    auto end()
    {
        return std::ranges::end(urange);
    }

    //!\copydoc end()
    auto end() const
    //!\cond
        requires seqan3::const_iterable_range<urng_t>
    //!\endcond
    {
        return std::ranges::end(urange);
    }

    //!\copydoc end()
    auto cend() const
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

    template <typename urng2_t>
    friend class window_iterator;

public:
    /*!\name Associated types
     * \{
     */
    //!\brief Type for distances between iterators.
    using difference_type = typename std::iter_difference_t<it_t>;
    //!\brief Value type of this iterator.
    using value_type = std::ranges::range_value_t<urng_t>;
    //!\brief The pointer type.
    using pointer = void;
    //!\brief Reference to `value_type`.
    using reference = value_type;
    /*!\brief Tag this class as a forward iterator if `urng_t` models std::ranges::forward range, otherwise as an
     *        input iterator.
     */
    using iterator_category = std::conditional_t<std::ranges::forward_iterator<it_t>,
                                                 std::forward_iterator_tag,
                                                 std::input_iterator_tag>;
    /*!\brief Tag this class as a forward iterator if `urng_t` models std::ranges::forward range, otherwise as an
     *        input iterator.
     */
    using iterator_concept = iterator_category;
    //!\}

    /*!\name Constructors, destructor and assignment
     * \{
     */
    window_iterator()                                    = default; //!< Defaulted.
    window_iterator(window_iterator const &)             = default; //!< Defaulted.
    window_iterator(window_iterator &&)                  = default; //!< Defaulted.
    window_iterator & operator=(window_iterator const &) = default; //!< Defaulted.
    window_iterator & operator=(window_iterator &&)      = default; //!< Defaulted.
    ~window_iterator()                                   = default; //!< Defaulted.

    //!\brief Allow iterator on a const range to be constructible from an iterator over a non-const range.
    template <typename urng2_t>
    //!\cond
        requires std::same_as<std::remove_const_t<urng_t>, urng2_t>
     //!\endcond
    window_iterator(window_iterator<urng2_t> it) :
         urange_end{std::move(it.urange_end)},
         minimiser_value{std::move(it.minimiser_value)},
         window_right{std::move(it.window_right)},
         window_values{std::move(it.window_values)}
    {}

    /*!\brief                              Construct from begin and end iterators of a given range over
    *                                      std::totally_ordered values, and the number of values per window.
    * /param[in] it_start                  Iterator pointing to the first position of the std::totally_ordered range.
    * /param[in] it_end                    Iterator pointing to the last position of the std::totally_ordered range.
    * /param[in] window_values_size        The number of values in one window.
    *
    * \details
    *
    * Looks at the number of values per window, returns the smallest as minimiser and shifts then by one to repeat
    * this action. If a minimiser in consecutive windows is the same, it is returned only once.
    *
    */
    window_iterator(it_t it_start, sentinel_t it_end, uint32_t window_values_size) :
    urange_end{it_end}, window_right{it_start}
    {
        if (window_values_size > std::ranges::distance(window_right, urange_end))
            window_values_size = std::ranges::distance(window_right, urange_end);
        if (window_right != urange_end)
            window_first(window_values_size);
    }
    //!\}

    //!\anchor window_iterator_comparison
    //!\name Comparison operators
    //!\{

    //!\brief Compare to the sentinel of the underlying range.
    friend bool operator==(window_iterator const & lhs, sentinel_t const & rhs)
    {
        return lhs.window_right == rhs;
    }

    //!\brief Compare to the sentinel of the underlying range.
    friend bool operator==(sentinel_t const & lhs, window_iterator const & rhs)
    {
        return lhs == rhs.window_right;
    }

    //!\brief Compare to another window_iterator.
    friend bool operator==(window_iterator const & lhs, window_iterator const & rhs)
    {
        return (lhs.window_right == rhs.window_right) &&
               (lhs.window_values.size() == rhs.window_values.size());
    }

    //!\brief Compare to the sentinel of the underlying range.
    friend bool operator!=(window_iterator const & lhs, sentinel_t const & rhs)
    {
        return !(lhs == rhs);
    }

    //!\brief Compare to the sentinel of the underlying range.
    friend bool operator!=(sentinel_t const & lhs, window_iterator const & rhs)
    {
        return !(lhs == rhs);
    }

    //!\brief Compare to another window_iterator.
    friend bool operator!=(window_iterator const & lhs, window_iterator const & rhs)
    {
        return !(lhs == rhs);
    }

    //!\brief Compare to another window_iterator.
    friend bool operator<(window_iterator const & lhs, window_iterator const & rhs)
    {
        return (lhs.window_right < rhs.window_right) && (lhs.window_values.size() < rhs.window_values.size());
    }

    //!\brief Compare to another window_iterator.
    friend bool operator>(window_iterator const & lhs, window_iterator const & rhs)
    {
        return (lhs.window_right > rhs.window_right) && (lhs.window_values.size() > rhs.window_values.size());
    }

    //!\brief Compare to another window_iterator.
    friend bool operator<=(window_iterator const & lhs, window_iterator const & rhs)
    {
        return (lhs.window_right <= rhs.window_right) && (lhs.window_values.size() <= rhs.window_values.size());
    }

    //!\brief Compare to another window_iterator.
    friend bool operator>=(window_iterator const & lhs, window_iterator const & rhs)
    {
        return (lhs.window_right >= rhs.window_right) && (lhs.window_values.size() >= rhs.window_values.size());
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

    //!\brief Iterator to the rightmost value of one window.
    it_t window_right;

    //!\brief Stored values per window. It is necessary to store them, because a shift can remove the current minimiser.
    std::deque<uint64_t> window_values;

    //!\brief Increments iterator by 1.
    void get_minimiser()
    {
        // Call next_minimiser until minimiser value changed or end of the underlying range is reached.
        while (!next_minimiser()) {}
    }

    //!\brief Calculates minimisers for the first window.
    void window_first(uint32_t window_values_size)
    {
        for (uint32_t i = 0; (i < window_values_size - 1) ; i++)
        {
            window_values.push_back(*window_right);
            std::ranges::advance(window_right,  1);
        }
        window_values.push_back(*window_right);
        minimiser_value = *(std::min_element(std::begin(window_values), std::end(window_values)));
    }

    //!\brief Calculates the next minimiser value.
    // For the following windows, we remove the first window value (is now not in window_values) and add the new
    // value that results from the window shifting.
    bool next_minimiser()
    {
        std::ranges::advance(window_right, 1);
        if (window_right == urange_end)
            return true;

        uint64_t new_value = *window_right;
        if (minimiser_value == *(std::begin(window_values)))
        {
            window_values.pop_front();
            window_values.push_back(new_value);
            minimiser_value = *(std::min_element(std::begin(window_values), std::end(window_values)));
            return true;
        }

        window_values.pop_front();
        window_values.push_back(new_value);

        if (new_value < minimiser_value)
        {
            minimiser_value = new_value;
            return true;
        }

        return false;
    }
};

//!\brief A deduction guide for the view class template.
template <std::ranges::viewable_range rng_t>
minimiser_view(rng_t &&, uint32_t const & window_values_size) -> minimiser_view<std::ranges::all_view<rng_t>>;


// ---------------------------------------------------------------------------------------------------------------------
// minimiser_fn (adaptor definition)
// ---------------------------------------------------------------------------------------------------------------------

//![adaptor_def]
//!\brief views::minimiser's range adaptor object type (non-closure).
struct minimiser_fn
{
    //!\brief Store the number of values in one window and return a range adaptor closure object.
    constexpr auto operator()(uint32_t const & window_values_size) const
    {
        return adaptor_from_functor{*this, window_values_size};
    }

    /*!\brief                           Call the view's constructor with two arguments: the underlying view and an
     *                                  integer indicating how many values one window contains.
     * \param[in] urange                The input range to process. Must model std::ranges::viewable_range and
     *                                  std::ranges::forward_range.
     * \param[in] window_values_size    The number of values in one window.
     * \returns                         A range of converted values.
     */
    template <std::ranges::range urng_t>
    //!\cond
        requires std::convertible_to<std::ranges::range_value_t<urng_t>, uint64_t>
    //!\endcond
    constexpr auto operator()(urng_t && urange, uint32_t const & window_values_size) const
    {
        static_assert(std::ranges::viewable_range<urng_t>,
            "The range parameter to views::minimiser cannot be a temporary of a non-view range.");
        static_assert(std::ranges::forward_range<urng_t>,
            "The range parameter to views::minimiser must model std::ranges::forward_range.");

        if (window_values_size == 1)
            throw std::invalid_argument{"The chosen window_valzes_size is not valid. "
                                        "Please choose a value greater than 1."};

        return minimiser_view{urange, window_values_size};
    }
};
//![adaptor_def]

} // namespace seqan3::detail

namespace seqan3::views
{

/*!\name General purpose views
 * \{
 */

/*!\brief                           Computes minimisers for a range of comparable values. A minimiser is the smallest
 *                                  value in a window.
 * \tparam urng_t                   The type of the range being processed. See below for requirements. [template
 *                                  parameter is omitted in pipe notation]
 * \param[in] urange                The range being processed. [parameter is omitted in pipe notation]
 * \param[in] window_values_size    The number of values in one window.
 * \returns                         A range of std::totally_ordered where each value is the minimal value for one window.
 *                                  See below for the properties of the returned range.
 * \ingroup views
 *
 * \details
 *
 * A minimiser is the smallest value in a window. For example for the following list of hash values
 * [28, 100, 9, 23, 4, 1, 72, 37, 8] and 4 as `window_values_size`, the minimiser values are [9,4,1].
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
