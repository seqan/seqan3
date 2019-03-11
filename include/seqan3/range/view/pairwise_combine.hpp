// -----------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2019, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2019, Knut Reinert & MPI für molekulare Genetik
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at: https://github.com/seqan/seqan3/blob/master/LICENSE
// -----------------------------------------------------------------------------------------------------

/*!\file
 * \brief Provides seqan3::view::pairwise_combine.
 * \author Rene Rahn <rene.rahn AT fu-berlin.de>
 */

#pragma once

#include <cmath>

#include <seqan3/range/concept.hpp>
#include <seqan3/range/view/detail.hpp>
#include <seqan3/std/ranges>
#include <seqan3/std/view/view_all.hpp>

namespace seqan3::detail
{
/*!\brief Generates all pairwise combinations of the elements in the underlying range.
 * \ingroup view
 * \tparam underlying_range_type The type of the underlying range; must model std::ranges::View,
 *                               std::ranges::ForwardRange and std::ranges::CommonRange.
 *
 * \details
 *
 * This view will provide a convenient way to iterate over all pairwise combinations of the elements of the underlying
 * range (in no defined order). A underlying range with `n` elements will therefore have `n choose 2` possible
 * combinations.
 *
 * ### Example
 *
 * \include test/snippet/range/view/pairwise_combine.cpp
 */
template <std::ranges::View underlying_range_type>
//!\cond
    requires  std::ranges::ForwardRange<underlying_range_type> && std::ranges::CommonRange<underlying_range_type>
// !\endcond
class pairwise_combine_view : public ranges::view_interface<pairwise_combine_view<underlying_range_type>>
{
private:

    //!\brief Alias type for the iterator over the underlying range.
    using underlying_iterator_type = std::ranges::iterator_t<underlying_range_type>;

    /*!\brief The internal iterator type.
     *
     * \details
     *
     * This iterator models the iterator category of the `underlying_iterator_type`. It maintains a pair of iterators on
     * the underlying range that move over all pairwise combinations of elements. The end is reached, when the second
     * iterator points to the end of the underlying range and the first iterator to the last element of the range.
     * Also note that this iterator does not model [Cpp17Iterator](https://en.cppreference.com/w/cpp/named_req/Iterator),
     * also known as `Cpp17Iterator` since it does not return a reference to the represented type but a prvalue.
     * Thus this iterator might not be usable with some legacy algorithms of the STL. But it is guaranteed to work with
     * the ranges algorithms.
     */
    class iterator_type
    {
    private:
        //!\brief Alias for the value type of the underlying iterator type.
        using underlying_val_t = typename std::iterator_traits<underlying_iterator_type>::value_type;
        //!\brief Alias for the reference type of the underlying iterator type.
        using underlying_ref_t = typename std::iterator_traits<underlying_iterator_type>::reference;
    public:

        /*!\name Associated types
         * \{
         */
        using difference_type   = std::ptrdiff_t;
        using value_type        = std::tuple<underlying_val_t, underlying_val_t>;
        using reference         = std::tuple<underlying_ref_t, underlying_ref_t>;
        using pointer           = void;
        using iterator_category = typename std::iterator_traits<underlying_iterator_type>::iterator_category;
        //!\}

        /*!\name Constructors, destructor and assignment
         * \{
         */
        constexpr iterator_type() = default;
        constexpr iterator_type(iterator_type const &) = default;
        constexpr iterator_type(iterator_type &&) = default;
        constexpr iterator_type & operator=(iterator_type const &) = default;
        constexpr iterator_type & operator=(iterator_type &&) = default;
        ~iterator_type() = default;

        /*!\brief Constructs the iterator from the current underlying iterator and the end iterator of the underlying
         *        range.
         * \param[in] iter     The iterator pointing to current element within the underlying range.
         * \param[in] begin_it The iterator pointing to begin of the underlying range. Only needed if
         *                     `underlying_iterator_type` models std::RandomAccessIterator.
         * \param[in] end_it   The iterator pointing to end of the underlying range.
         *
         * \details
         *
         * Constructs the iterator by caching the end of the underlying range and setting the first iterator to the
         * given position and the second to the first incremented by one.
         */
        constexpr iterator_type(underlying_iterator_type iter,
                                underlying_iterator_type begin_it,
                                underlying_iterator_type end_it) noexcept :
            first_it{iter},
            second_it{++iter},
            begin_it{begin_it},
            end_it{end_it}
        {}
        //!\}

        /*!\name Accessors
         * \{
         */
         //!\brief Accesses the pointed-to element
        constexpr reference operator*() const
            noexcept(noexcept(*std::declval<underlying_iterator_type>()))
        {
            return {*first_it, *second_it};
        }

        /*!\brief Access the element at the given index
         * \param[in] index The index of the element to be returned.
         */
        constexpr reference operator[](size_t const index)
            noexcept(noexcept(std::declval<iterator_type &>().from_index(1)))
        //!\cond
            requires std::RandomAccessIterator<underlying_iterator_type>
        //!\endcond
        {
            from_index(index);
            return **this;
        }
        //!\}

        /*!\name Arithmetic operators
         * \{
         */
        //!\brief Pre-increment operator.
        constexpr iterator_type & operator++(/*pre-increment*/)
            noexcept(noexcept(++std::declval<underlying_iterator_type &>()))
        {
            if (++second_it == end_it)
            {
                ++first_it;
                second_it = first_it;
                ++second_it;
            }
            return *this;
        }

        //!\brief Post-increment operator.
        constexpr iterator_type operator++(int /*post-increment*/)
            noexcept(noexcept(std::declval<underlying_iterator_type &>()++))
        {
            iterator_type tmp{*this};
            ++*this;
            return tmp;
        }

        //!\brief Pre-decrement operator; `underlying_iterator_type` must model std::BidirectionalIterator.
        constexpr iterator_type & operator--(/*pre-decrement*/)
            noexcept(noexcept(--std::declval<underlying_iterator_type &>()))
        //!\cond
            requires std::BidirectionalIterator<underlying_iterator_type>
        //!\endcond
        {
            if (--second_it == first_it)
            {
                --first_it;
                second_it = end_it;
                --second_it;
            }
            return *this;
        }

        //!\brief Post-decrement operator; `underlying_iterator_type` must model std::BidirectionalIterator.
        constexpr iterator_type operator--(int /*post-decrement*/)
            noexcept(noexcept(std::declval<underlying_iterator_type &>()--))
        //!\cond
            requires std::BidirectionalIterator<underlying_iterator_type>
        //!\endcond
        {
            iterator_type tmp{*this};
            --*this;
            return tmp;
        }

        //!\brief Advances the iterator by the given offset; `underlying_iterator_type` must model
        //!\      std::RandomAccessIterator.
        constexpr iterator_type & operator+=(difference_type const offset)
            noexcept(noexcept(std::declval<iterator_type &>().from_index(1)))
        //!\cond
            requires std::RandomAccessIterator<underlying_iterator_type>
        //!\endcond
        {
            from_index(to_index() + offset);
            return *this;
        }

        //!\brief Advances the iterator by the given offset; `underlying_iterator_type` must model
        //!\      std::RandomAccessIterator.
        constexpr iterator_type operator+(difference_type const offset)
            noexcept(noexcept(std::declval<iterator_type &>() += 1))
        //!\cond
            requires std::RandomAccessIterator<underlying_iterator_type>
        //!\endcond
        {
            iterator_type tmp{*this};
            return (tmp += offset);
        }

        //!\brief Advances the iterator by the given offset; `underlying_iterator_type` must model
        //!\      std::RandomAccessIterator.
        constexpr friend iterator_type operator+(difference_type const offset, iterator_type iter)
            noexcept(noexcept(std::declval<iterator_type &>().from_index(1)))
        //!\cond
            requires std::RandomAccessIterator<underlying_iterator_type>
        //!\endcond
        {
            iter.from_index(iter.to_index() + offset);
            return iter;
        }

        //!\brief Decrements the iterator by the given offset; `underlying_iterator_type` must model
        //!\      std::RandomAccessIterator.
        constexpr iterator_type & operator-=(difference_type const offset)
            noexcept(noexcept(std::declval<iterator_type &>().from_index(1)))
        //!\cond
            requires std::RandomAccessIterator<underlying_iterator_type>
        //!\endcond
        {
            from_index(to_index() - offset);
            return *this;
        }

        //!\brief Decrements the iterator by the given offset; `underlying_iterator_type` must model
        //!\      std::RandomAccessIterator.
        constexpr iterator_type operator-(difference_type const offset)
            noexcept(noexcept(std::declval<iterator_type &>() -= 1))
        //!\cond
            requires std::RandomAccessIterator<underlying_iterator_type>
        //!\endcond
        {
            iterator_type tmp{*this};
            return (tmp -= offset);
        }

        //!\brief Computes the distance between two iterators; `underlying_iterator_type` must model
        //!\      std::RandomAccessIterator.
        constexpr friend difference_type operator-(iterator_type const & lhs, iterator_type const & rhs)
            noexcept(noexcept(std::declval<iterator_type &>().to_index()))
        //!\cond
            requires std::RandomAccessIterator<underlying_iterator_type>
        //!\endcond
        {
            return static_cast<difference_type>(lhs.to_index() - rhs.to_index());
        }
        //!\}

        /*!\name Comparison operators
         * \brief These operators are available if the `underlying_iterator_type` models std::EqualityComparable or
         *        std::StrictTotallyOrdered respectively.
         * \{
         */
        constexpr friend bool operator==(iterator_type const & lhs, iterator_type const & rhs)
            noexcept(noexcept(std::declval<underlying_iterator_type &>() == std::declval<underlying_iterator_type &>()))
        //!\cond
            requires std::EqualityComparable<underlying_iterator_type>
        //!\endcond
        {
            return std::tie(lhs.first_it, lhs.second_it) == std::tie(rhs.first_it, rhs.second_it);
        }

        constexpr friend bool operator!=(iterator_type const & lhs, iterator_type const & rhs)
            noexcept(noexcept(std::declval<underlying_iterator_type &>() != std::declval<underlying_iterator_type &>()))
        //!\cond
            requires std::EqualityComparable<underlying_iterator_type>
        //!\endcond
        {
            return !(lhs == rhs);
        }

        constexpr friend bool operator<(iterator_type const & lhs, iterator_type const & rhs)
            noexcept(noexcept(std::declval<underlying_iterator_type &>() < std::declval<underlying_iterator_type &>()))
        //!\cond
            requires std::StrictTotallyOrdered<underlying_iterator_type>
        //!\endcond
        {
            return std::tie(lhs.first_it, lhs.second_it) < std::tie(rhs.first_it, rhs.second_it);
        }

        constexpr friend bool operator>(iterator_type const & lhs, iterator_type const & rhs)
            noexcept(noexcept(std::declval<underlying_iterator_type &>() > std::declval<underlying_iterator_type &>()))
        //!\cond
            requires std::StrictTotallyOrdered<underlying_iterator_type>
        //!\endcond
        {
            return std::tie(lhs.first_it, lhs.second_it) > std::tie(rhs.first_it, rhs.second_it);
        }

        constexpr friend bool operator<=(iterator_type const & lhs, iterator_type const & rhs)
            noexcept(noexcept(std::declval<underlying_iterator_type &>() <= std::declval<underlying_iterator_type &>()))
        //!\cond
            requires std::StrictTotallyOrdered<underlying_iterator_type>
        //!\endcond
        {
            return std::tie(lhs.first_it, lhs.second_it) <= std::tie(rhs.first_it, rhs.second_it);
        }

        constexpr friend bool operator>=(iterator_type const & lhs, iterator_type const & rhs)
            noexcept(noexcept(std::declval<underlying_iterator_type &>() >= std::declval<underlying_iterator_type &>()))
        //!\cond
            requires std::StrictTotallyOrdered<underlying_iterator_type>
        //!\endcond
        {
            return std::tie(lhs.first_it, lhs.second_it) >= std::tie(rhs.first_it, rhs.second_it);
        }
        //!\}

    private:

        /*!\brief Returns the index for the current iterator position.
         *
         * \details
         *
         * The pairwise combination can also be seen as a triangular matrix, where given a size n the corresponding
         * matrix has values in [0, 1], [0, 2], ... , [0, n - 1], [1, 2], ..., [1, n - 1], ... [n - 2, n - 1] and all
         * other entries are empty. Using this scheme one can use the properties of triangular numbers, where the
         * diagonal index of this matrix can be computed using triangular roots
         * (see https://stackoverflow.com/questions/27086195/linear-index-upper-triangular-matrix).
         * Given these properties, one can compute the matrix index (i, j) from the linearised matrix index and vice
         * versa.
         */
        constexpr size_t to_index() const
            noexcept(noexcept(std::declval<underlying_iterator_type &>() - std::declval<underlying_iterator_type &>()))
        //!\cond
            requires std::RandomAccessIterator<underlying_iterator_type>
        //!\endcond
        {
            size_t src_size = end_it - begin_it;
            size_t index_i = first_it - begin_it;
            size_t index_j = second_it - begin_it;
            return (src_size * (src_size - 1)/2) - (src_size - index_i) * ((src_size - index_i) - 1)/2 +
                   index_j - index_i - 1;
        }

        /*!\brief Sets the iterator to the given index.
         * \param[in] index The index to set the iterator to.
         * \copydetails to_index()
         */
        constexpr void from_index(size_t const index)
            noexcept(noexcept(std::declval<underlying_iterator_type &>() - std::declval<underlying_iterator_type &>()) &&
                     noexcept(std::declval<underlying_iterator_type &>() + 1))
        //!\cond
            requires std::RandomAccessIterator<underlying_iterator_type>
        //!\endcond
        {
            size_t src_size = end_it - begin_it;
            size_t index_i = src_size - 2 -
                             std::floor(std::sqrt(-8 * index + 4 * src_size * (src_size - 1) - 7)/2.0 - 0.5);
            size_t index_j = index + index_i + 1 - src_size * (src_size - 1)/2 + (src_size - index_i) *
                             ((src_size - index_i) - 1)/2;
            first_it = begin_it + index_i;
            second_it = begin_it + index_j;
        }

        //!\brief The iterator pointing to the first element of the pairwise combination.
        underlying_iterator_type first_it{};
        //!\brief The iterator pointing to the second element of the pairwise combination.
        underlying_iterator_type second_it{};
        //!\brief The begin of the underlying range.
        underlying_iterator_type begin_it{};
        //!\brief The end of the underlying range.
        underlying_iterator_type end_it{};
    };

public:

    /*!\name Constructors, destructor and assignment
     * \{
     */
    constexpr pairwise_combine_view() = default;
    constexpr pairwise_combine_view(pairwise_combine_view const &) = default;
    constexpr pairwise_combine_view(pairwise_combine_view &&) = default;
    constexpr pairwise_combine_view & operator=(pairwise_combine_view const &) = default;
    constexpr pairwise_combine_view & operator=(pairwise_combine_view &&) = default;
    ~pairwise_combine_view() = default;

    /*!\brief Constructs from a view.
     * \param[in] range The underlying range to be wrapped. Of type `underlying_range_type`.
     *
     * \details
     *
     * During construction the iterator pointing to the last element of the view is cached (not the end of the range).
     * This optimises the call to end if the underlying range models only ForwardRange. Otherwise the call to end will be
     * linear in the number of elements of the underlying range.
     *
     * \attention This view cannot be chained immediately after an infinite range, because upon construction it will
     *            take forever to reach the last element of the view.
     *
     * ### Complexity
     *
     * Constant if `underlying_range_type` models std::ranges::BidirectionalRange, otherwise linear.
     */
    constexpr pairwise_combine_view(underlying_range_type range) : src_range{std::move(range)}
    {
        // Check if range is empty.
        if (std::ranges::empty(src_range))
        {
            back_iterator = std::ranges::end(src_range);
        }
        else
        {
            if constexpr (std::ranges::BidirectionalRange<underlying_range_type>)
            { // Simply take one before the end. We can do this as we require underlying_range_type to be a common range.
                back_iterator = std::ranges::prev(std::ranges::end(src_range));
            }
            else
            { // For all other cases we need to set the back_iterator in linear time to the correct position.
                back_iterator = std::ranges::begin(src_range);
                if constexpr (std::ranges::SizedRange<underlying_range_type>)
                {
                    std::ranges::advance(back_iterator, std::ranges::size(src_range) - 1);
                }
                else // We don't have the size, so we need to increment until one before the end in a linear pass.
                {
                    auto tmp_it = back_iterator;
                    do
                    {
                        back_iterator = tmp_it;
                    } while (++tmp_it != std::ranges::end(src_range));
                }
            }
        }
    }

    /*!\brief Constructs from a view.
     * \tparam    other_range_t  The type of the range to be wrapped with seqan3::detail::pairwise_combine_view;
     *                           must model std::ranges::ViewableRange and underlying_range_type must be constructible
     *                           with other_range wrapped in view::all.
     * \param[in] range          The underlying range to be wrapped.
     *
     * \details
     *
     * During construction the iterator pointing to the last element of the view is cached (not the end of the range).
     * This optimises the call to end if the underlying range models only ForwardRange. Otherwise the call to end will be
     * linear in the number of elements of the underlying range.
     *
     * \attention This view cannot be chained immediately after an infinite range, because upon construction it will
     *            take forever to reach the last element of the view.
     *
     * ### Complexity
     *
     * Constant if `other_range_t` models std::ranges::BidirectionalRange, otherwise linear.
     */
    template <typename other_range_t>
    //!\cond
        requires !std::Same<remove_cvref_t<other_range_t>, pairwise_combine_view> &&
                 std::ranges::ViewableRange<other_range_t> &&  // Must come after self type check to avoid conflicts with the move constructor.
                 std::Constructible<underlying_range_type, ranges::ref_view<std::remove_reference_t<other_range_t>>>
            //TODO: Investigate: the following expression is equivalent to the one above but raises a weird assertion in
            //      the ranges adaptor suggesting that the pairwise_combine_view is not a ViewableRange.
            //      std::Constructible<underlying_range_type, decltype(view::all(std::declval<other_range_t &&>()))>
    //!\endcond
    constexpr pairwise_combine_view(other_range_t && range) :
        pairwise_combine_view{view::all(std::forward<other_range_t>(range))}
    {}

    /*!\name Iterators
     * \{
     */
    /*!\brief Returns an iterator to the first element of the range.
     * \returns Iterator to the first element.
     *
     * If the range is empty, the returned iterator will be equal to end().
     *
     * ### Complexity
     *
     * Constant.
     *
     * ### Exceptions
     *
     * No-throw guarantee.
     */
    constexpr iterator_type begin() noexcept
    {
        return iterator_type{std::ranges::begin(src_range), std::ranges::begin(src_range), std::ranges::end(src_range)};
    }

    //!\copydoc begin()
    constexpr iterator_type begin() const noexcept
    {
        return iterator_type{std::ranges::begin(src_range), std::ranges::begin(src_range), std::ranges::end(src_range)};
    }

    //!\copydoc begin()
    constexpr iterator_type cbegin() const noexcept
    {
        return begin();
    }

    /*!\brief Returns an iterator to the element following the last element of the range.
     * \returns Iterator to the end.
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
    constexpr iterator_type end() noexcept
    {
        return iterator_type{back_iterator, std::ranges::begin(src_range), std::ranges::end(src_range)};
    }

    //!\copydoc end()
    constexpr iterator_type end() const noexcept
    {
        return iterator_type{back_iterator, std::ranges::begin(src_range), std::ranges::end(src_range)};
    }

    //!\copydoc end()
    constexpr iterator_type cend() const noexcept
    {
        return end();
    }
    //!\}

private:

    //!\brief The underling range.
    underlying_range_type src_range{};
    //!\brief The cached iterator pointing to the last element of the underlying range.
    std::ranges::iterator_t<underlying_range_type> back_iterator{};
};

/*!\name Type deduction guides
 * \{
 */

//!\brief Deduces the correct template type from a view.
template <std::ranges::View other_range_t>
pairwise_combine_view(other_range_t range) -> pairwise_combine_view<other_range_t>;

//!\brief Deduces the correct template type from a non-view lvalue range by wrapping the range in seqan3::view::all.
template <std::ranges::ViewableRange other_range_t>
pairwise_combine_view(other_range_t && range) ->
    pairwise_combine_view<decltype(view::all(std::declval<other_range_t &&>()))>;
//!\}

} // namespace seqan3::detail

namespace seqan3::view
{
/*!\name General purpose views
 * \{
 */

/*!\brief             A view adapter that generates all pairwise combinations of the elements of the underlying range.
 * \tparam urng_t     The type of the range being processed. See below for requirements.
 * \param[in] urange  The range being processed.
 * \returns           A view over all pairwise combinations of the elements of the underlying range.
 *                    See below for the properties of the returned range.
 * \ingroup view
 *
 * \details
 *
 * This view generates two-element tuples representing all unique combinations of the elements of the underlying range
 * (the order of the elements does not matter). If the underlying range has less than two elements the returned range is
 * empty, otherwise the size of the returned range corresponds to the binomial coefficient `n choose 2`, where `n` is
 * the size of the underlying range. The reference type of this range is a tuple over the reference type of the
 * underlying range. This range is const-iterable as long as the underlying range is iterable (const-ness of this range
 * will not be propagated to the underlying range).
 * In order to receive the end iterator in constant time an iterator pointing to the last element of the underlying
 * range will be cached upon construction of this view. This construction takes linear time for underlying ranges that
 * do not model std::ranges::BidirectionalRange.
 *
 * ### Iterator
 *
 * The returned iterator from begin does not model std::LegacyIterator, also known as `Cpp17Iterator` since it
 * does not return a reference to the represented type but a prvalue. Thus this iterator might not be usable within
 * some legacy algorithms of the STL. But it is guaranteed to work with the ranges algorithms.
 *
 * ### View properties
 *
 * | range concepts and reference_t  | `urng_t` (underlying range type)      | `rrng_t` (returned range type)                                       |
 * |---------------------------------|:-------------------------------------:|:--------------------------------------------------------------------:|
 * | std::ranges::InputRange         | *required*                            | *preserved*                                                          |
 * | std::ranges::ForwardRange       | *required*                            | *preserved*                                                          |
 * | std::ranges::BidirectionalRange |                                       | *preserved*                                                          |
 * | std::ranges::RandomAccessRange  |                                       | *preserved*                                                          |
 * | std::ranges::ContiguousRange    |                                       | *lost*                                                               |
 * |                                 |                                       |                                                                      |
 * | std::ranges::ViewableRange      | *required*                            | *guaranteed*                                                         |
 * | std::ranges::View               |                                       | *guaranteed*                                                         |
 * | std::ranges::SizedRange         |                                       | *preserved*                                                          |
 * | std::ranges::CommonRange        | *required*                            | *guaranteed*                                                         |
 * | std::ranges::OutputRange        |                                       | *lost*                                                               |
 * | seqan3::const_iterable_concept  |                                       | *preserved*                                                          |
 * |                                 |                                       |                                                                      |
 * | seqan3::reference_t             |                                       | std::tuple<seqan3::reference_t<urng_t>, seqan3::reference_t<urng_t>> |
 *
 * See the \link view view submodule documentation \endlink for detailed descriptions of the view properties.
 *
 * ### Thread safety
 *
 * Concurrent access to this view, e.g. while iterating over it, is thread-safe and must not be protected externally.
 *
 * ### Example
 *
 * \include test/snippet/range/view/pairwise_combine.cpp
 *
 * \attention This view cannot be chained immediately with an infinite range, because upon construction it will
 *            take forever to reach the last element of the view.
 *
 * \hideinitializer
 */
inline constexpr auto pairwise_combine = detail::generic_pipable_view_adaptor<detail::pairwise_combine_view>{};

//!\}
} // namespace seqan3::view
