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
/*!\brief Generates all pairwise combinations of the elements in the source range.
 * \ingroup view
 * \tparam source_range_type The type of the source range; must model std::ranges::View, std::ranges::ForwardRange
 *
 * \details
 *
 * This view will provide a convenient way to iterate over all pairwise combinations of a given source range
 * while the order of selection does not matter. A source range with `n` elements will therefore have `n choose 2`
 * possible combinations.
 */
template <std::ranges::View source_range_type>
//!\cond
    requires  std::ranges::ForwardRange<source_range_type>
// !\endcond
class pairwise_combine_view : public ranges::view_interface<pairwise_combine_view<source_range_type>>
{
private:

    //!\brief Alias type for the iterator over the source range.
    using source_iterator_type = std::ranges::iterator_t<source_range_type>;

    /*!\brief The internal iterator type.
     *
     * \details
     *
     * This iterator models the iterator category of the `source_iterator_type`. It maintains a pair of iterators on
     * the source range that move over all pairwise combinations of elements. The end is reached, when the second
     * iterator points to the end of the source range and the first iterator to the last element of the range.
     */
    class iterator_type
    {
    private:
        //!\brief Alias for the value type of the source iterator type.
        using source_val_t = typename std::iterator_traits<source_iterator_type>::value_type;
        //!\brief Alias for the reference type of the source iterator type.
        using source_ref_t = typename std::iterator_traits<source_iterator_type>::reference;
    public:

        /*!\name Associated types
         * \{
         */
        using difference_type   = std::ptrdiff_t;
        using value_type        = std::tuple<source_val_t, source_val_t>;
        using reference         = std::tuple<source_ref_t, source_ref_t>;
        using pointer           = void;
        using iterator_category = typename std::iterator_traits<source_iterator_type>::iterator_category;
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

        /*!\brief Constructs the iterator from the current source iterator and the end iterator of the source range.
         * \param[in] iter     The iterator pointing to current element within the source range.
         * \param[in] begin_it The iterator pointing to begin of the source range. Only needed if
         *                     `source_iterator_type` models std::RandomAccessIterator.
         * \param[in] end_it   The iterator pointing to end of the source range.
         *
         * \details
         *
         * Constructs the iterator by caching the end of the source range and setting the first iterator to the
         * given position and the second to the first incremented by one.
         */
        constexpr iterator_type(source_iterator_type iter,
                                source_iterator_type begin_it,
                                source_iterator_type end_it) noexcept :
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
            noexcept(noexcept(*std::declval<source_iterator_type>()))
        {
            return {*first_it, *second_it};
        }

        /*!\brief Access the element at the given index
         * \param[in] index The index of the element to be returned.
         */
        constexpr reference operator[](size_t const index)
            noexcept(noexcept(std::declval<iterator_type &>().from_index(1)))
        //!\cond
            requires std::RandomAccessIterator<source_iterator_type>
        //!\endcond
        {
            from_index(index);
            return this->operator*();
        }
        //!\}

        /*!\name Arithmetic operators
         * \{
         */
        //!\brief Pre-increment operator.
        constexpr iterator_type & operator++(/*pre-increment*/)
            noexcept(noexcept(++std::declval<source_iterator_type &>()))
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
            noexcept(noexcept(std::declval<source_iterator_type &>()++))
        {
            iterator_type tmp{*this};
            ++*this;
            return tmp;
        }

        //!\brief Pre-decrement operator; `source_iterator_type` must model std::BidirectionalIterator.
        constexpr iterator_type & operator--(/*pre-decrement*/)
            noexcept(noexcept(--std::declval<source_iterator_type &>()))
        //!\cond
            requires std::BidirectionalIterator<source_iterator_type>
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

        //!\brief Post-decrement operator; `source_iterator_type` must model std::BidirectionalIterator.
        constexpr iterator_type operator--(int /*post-decrement*/)
            noexcept(noexcept(std::declval<source_iterator_type &>()--))
        //!\cond
            requires std::BidirectionalIterator<source_iterator_type>
        //!\endcond
        {
            iterator_type tmp{*this};
            --*this;
            return tmp;
        }

        //!\brief Advances the iterator by the given offset; `source_iterator_type` must model std::RandomAccessIterator.
        constexpr iterator_type & operator+=(difference_type const offset)
            noexcept(noexcept(std::declval<iterator_type &>().from_index(1)))
        //!\cond
            requires std::RandomAccessIterator<source_iterator_type>
        //!\endcond
        {
            from_index(to_index() + offset);
            return *this;
        }

        //!\brief Advances the iterator by the given offset; `source_iterator_type` must model std::RandomAccessIterator.
        constexpr iterator_type operator+(difference_type const offset)
            noexcept(noexcept(std::declval<iterator_type &>() += 1))
        //!\cond
            requires std::RandomAccessIterator<source_iterator_type>
        //!\endcond
        {
            iterator_type tmp{*this};
            return (tmp += offset);
        }

        //!\brief Advances the iterator by the given offset; `source_iterator_type` must model std::RandomAccessIterator.
        constexpr friend iterator_type operator+(difference_type const offset, iterator_type iter)
            noexcept(noexcept(std::declval<iterator_type &>().from_index(1)))
        //!\cond
            requires std::RandomAccessIterator<source_iterator_type>
        //!\endcond
        {
            iter.from_index(iter.to_index() + offset);
            return iter;
        }

        //!\brief Decrements the iterator by the given offset; `source_iterator_type` must model std::RandomAccessIterator.
        constexpr iterator_type & operator-=(difference_type const offset)
            noexcept(noexcept(std::declval<iterator_type &>().from_index(1)))
        //!\cond
            requires std::RandomAccessIterator<source_iterator_type>
        //!\endcond
        {
            from_index(to_index() - offset);
            return *this;
        }

        //!\brief Decrements the iterator by the given offset; `source_iterator_type` must model std::RandomAccessIterator.
        constexpr iterator_type operator-(difference_type const offset)
            noexcept(noexcept(std::declval<iterator_type &>() -= 1))
        //!\cond
            requires std::RandomAccessIterator<source_iterator_type>
        //!\endcond
        {
            iterator_type tmp{*this};
            return (tmp -= offset);
        }

        //!\brief Computes the distance between two iterators; `source_iterator_type` must model std::RandomAccessIterator.
        constexpr friend difference_type operator-(iterator_type const & lhs, iterator_type const & rhs)
            noexcept(noexcept(std::declval<iterator_type &>().to_index()))
        //!\cond
            requires std::RandomAccessIterator<source_iterator_type>
        //!\endcond
        {
            return static_cast<difference_type>(lhs.to_index() - rhs.to_index());
        }
        //!\}

        /*!\name Comparison operators
         * \brief These operators are available if the `source_iterator_type` models std::EqualityComparable or
         *        std::StrictTotallyOrdered respectively.
         * \{
         */
        constexpr friend bool operator==(iterator_type const & lhs, iterator_type const & rhs)
            noexcept(noexcept(std::declval<source_iterator_type &>() == std::declval<source_iterator_type &>()))
        //!\cond
            requires std::EqualityComparable<source_iterator_type>
        //!\endcond
        {
            return std::tie(lhs.first_it, lhs.second_it) == std::tie(rhs.first_it, rhs.second_it);
        }

        constexpr friend bool operator!=(iterator_type const & lhs, iterator_type const & rhs)
            noexcept(noexcept(std::declval<source_iterator_type &>() != std::declval<source_iterator_type &>()))
        //!\cond
            requires std::EqualityComparable<source_iterator_type>
        //!\endcond
        {
            return !(lhs == rhs);
        }

        constexpr friend bool operator<(iterator_type const & lhs, iterator_type const & rhs)
            noexcept(noexcept(std::declval<source_iterator_type &>() < std::declval<source_iterator_type &>()))
        //!\cond
            requires std::StrictTotallyOrdered<source_iterator_type>
        //!\endcond
        {
            return std::tie(lhs.first_it, lhs.second_it) < std::tie(rhs.first_it, rhs.second_it);
        }

        constexpr friend bool operator>(iterator_type const & lhs, iterator_type const & rhs)
            noexcept(noexcept(std::declval<source_iterator_type &>() > std::declval<source_iterator_type &>()))
        //!\cond
            requires std::StrictTotallyOrdered<source_iterator_type>
        //!\endcond
        {
            return std::tie(lhs.first_it, lhs.second_it) > std::tie(rhs.first_it, rhs.second_it);
        }

        constexpr friend bool operator<=(iterator_type const & lhs, iterator_type const & rhs)
            noexcept(noexcept(std::declval<source_iterator_type &>() <= std::declval<source_iterator_type &>()))
        //!\cond
            requires std::StrictTotallyOrdered<source_iterator_type>
        //!\endcond
        {
            return std::tie(lhs.first_it, lhs.second_it) <= std::tie(rhs.first_it, rhs.second_it);
        }

        constexpr friend bool operator>=(iterator_type const & lhs, iterator_type const & rhs)
            noexcept(noexcept(std::declval<source_iterator_type &>() >= std::declval<source_iterator_type &>()))
        //!\cond
            requires std::StrictTotallyOrdered<source_iterator_type>
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
            noexcept(noexcept(std::ranges::distance(std::declval<source_iterator_type &>(),
                                                    std::declval<source_iterator_type &>())))
        //!\cond
            requires std::RandomAccessIterator<source_iterator_type>
        //!\endcond
        {
            size_t src_size = std::ranges::distance(begin_it, end_it);
            size_t index_i = std::ranges::distance(begin_it, first_it);
            size_t index_j = std::ranges::distance(begin_it, second_it);
            return (src_size * (src_size - 1)/2) - (src_size - index_i) * ((src_size - index_i) - 1)/2 +
                   index_j - index_i - 1;
        }

        /*!\brief Sets the iterator to the given index.
         * \param[in] index The index to set the iterator to.
         * \copydetails to_index()
         */
        constexpr void from_index(size_t const index)
            noexcept(noexcept(std::declval<source_iterator_type &>() - std::declval<source_iterator_type &>()) &&
                     noexcept(std::declval<source_iterator_type &>() + 1))
        //!\cond
            requires std::RandomAccessIterator<source_iterator_type>
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
        source_iterator_type first_it{};
        //!\brief The iterator pointing to the second element of the pairwise combination.
        source_iterator_type second_it{};
        //!\brief The begin of the source range.
        source_iterator_type begin_it{};
        //!\brief The end of the source range.
        source_iterator_type end_it{};
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
     * \param[in] range A rvalue reference to the source view.
     *
     * \details
     *
     * During construction the iterator pointing to the last element of the view is cached. This optimises
     * the call to end if the source range models only ForwardRange. Otherwise the call to end will be
     * linear in the number of elements of the source range.
     */
    explicit constexpr pairwise_combine_view(source_range_type && range) : src_range{std::move(range)}
    {
        // Check if range is empty.
        if (std::ranges::begin(src_range) == std::ranges::end(src_range))
        {
            back_iterator = std::ranges::end(src_range);
        }
        else
        {
            if constexpr (std::ranges::BidirectionalRange<source_range_type> &&
                          std::ranges::CommonRange<source_range_type>)
            {
                back_iterator = std::ranges::prev(std::ranges::end(src_range));
            }
            else
            {  // Cannot decrement so we invoke linear pass until we get to the last element.
                size_t size = std::distance(src_range);
                back_iterator = std::ranges::begin(src_range);
                std::advance(back_iterator, size - 1);
            }
        }
    }

    /*!\brief Constructs from a range.
     * \param[in] range A lvalue reference to the source range.
     *
     * \details
     *
     * Wraps the passed range in a ref_view in order to store lvalue references to ranges (non-views) as a view
     * and delegate to the view constructor above.
     */
    template <typename other_range_type>
    explicit constexpr pairwise_combine_view(other_range_type & range) : pairwise_combine_view{view::all(range)}
    {}
    //!\}

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
        return iterator_type{std::ranges::begin(src_range), std::ranges::begin(src_range), std::ranges::end(src_range)};
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
        return iterator_type{back_iterator, std::ranges::begin(src_range), std::ranges::end(src_range)};
    }
    //!\}

    /*!\name Capacity
     * \{
     */

    /*!\brief Returns the size of the range.
     *
     * \details
     *
     * The size can only be computed when the entire range models std::ranges::RandomAccessRange.
     *
     * ### Complexity
     *
     * Constant.
     *
     * ### Exceptions
     *
     * No-throw guarantee.
     */
    constexpr size_t size() const noexcept
    //!\cond
        requires std::RandomAccessIterator<iterator_type>
    //!\endcond
    {
        return static_cast<size_t>(cend() - cbegin());
    }
    //!\}

private:

    //!\brief The underling range.
    source_range_type src_range{};
    //!\brief The cached iterator pointing to the last element of the source range.
    std::ranges::iterator_t<source_range_type> back_iterator{};
};

/*!\name Type deduction guides
 * \relates seqan3::detail::pairwise_combine_view
 * \{
 */

//!\brief Deduces the source range type from the passed view.
template <std::ranges::View other_type>
pairwise_combine_view(other_type view) -> pairwise_combine_view<std::remove_reference_t<other_type>>;

//!\brief Deduces the source range type from the passed range.
template <std::ranges::Range range_type>
    requires !std::ranges::View<range_type>
pairwise_combine_view(range_type & range) -> pairwise_combine_view<decltype(std::declval<range_type &>() | view::all)>;
//!\}

} // namespace seqan3::detail

namespace seqan3::view
{
/*!\name General purpose views
 * \{
 */

/*!\brief             A view adapter that generates all pairwise combinations of the elements of the source range.
 * \tparam urng_t     The type of the range being processed. See below for requirements.
 * \param[in] urange  The range being processed.
 * \returns           A range over all pairwise combinations. See below for the properties of the returned range.
 * \ingroup view
 *
 * \details
 *
 * This view generates two-element tuples representing all possible combinations of the elements of the source range
 * while ignoring the order of the elements. If the source range has less than two elements the returned range is empty,
 * otherwise the size of the returned view corresponds to the binomial coefficient `n choose 2`, where `n` is the
 * size of the source range. The reference type of this range is a tuple over the reference type of the source range.
 * Constness will not be propagated to the reference types of the source range.
 *
 * ### View properties
 *
 * | range concepts and reference_t  | `urng_t` (underlying range type)      | `rrng_t` (returned range type)                                       |
 * |---------------------------------|:-------------------------------------:|:--------------------------------------------------------------------:|
 * | std::ranges::InputRange         |                                       | *undefined*                                                          |
 * | std::ranges::ForwardRange       | *required*                            | *preserved*                                                          |
 * | std::ranges::BidirectionalRange |                                       | *preserved*                                                          |
 * | std::ranges::RandomAccessRange  |                                       | *preserved*                                                          |
 * | std::ranges::ContiguousRange    |                                       | *lost*                                                               |
 * |                                 |                                       |                                                                      |
 * | std::ranges::ViewableRange      | *required*                            | *guaranteed*                                                         |
 * | std::ranges::View               |                                       | *guaranteed*                                                         |
 * | std::ranges::SizedRange         |                                       | *guaranteed* iff urng_t models std::RandomAccessRange                |
 * | std::ranges::CommonRange        |                                       | *guaranteed*                                                         |
 * | std::ranges::OutputRange        |                                       | *lost*                                                               |
 * | seqan3::const_iterable_concept  |                                       | *preserved*                                                          |
 * |                                 |                                       |                                                                      |
 * | seqan3::reference_t             |                                       | std::tuple<seqan3::reference_t<urng_t>, seqan3::reference_t<urng_t>> |
 *
 * See the \link view view submodule documentation \endlink for detailed descriptions of the view properties.
 *
 * ### Thread safety
 *
 * Concurrent access to this view, e.g. while iterating over it, is not thread-safe and must be protected externally.
 *
 * \hideinitializer
 */
inline constexpr auto pairwise_combine = detail::generic_pipable_view_adaptor<detail::pairwise_combine_view>{};

//!\}
} // namespace seqan3::view
