// -----------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2020, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2020, Knut Reinert & MPI für molekulare Genetik
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at: https://github.com/seqan/seqan3/blob/master/LICENSE.md
// -----------------------------------------------------------------------------------------------------

/*!\file
 * \brief Provides seqan3::views::pairwise_combine.
 * \author Rene Rahn <rene.rahn AT fu-berlin.de>
 */

#pragma once

#include <cmath>

#include <seqan3/core/common_tuple.hpp>
#include <seqan3/range/concept.hpp>
#include <seqan3/range/views/detail.hpp>
#include <seqan3/std/ranges>

namespace seqan3::detail
{
/*!\brief Generates all pairwise combinations of the elements in the underlying range.
 * \ingroup views
 * \tparam underlying_range_type The type of the underlying range; must model std::ranges::view,
 *                               std::ranges::forward_range and std::ranges::common_range.
 *
 * \details
 *
 * This view will provide a convenient way to iterate over all pairwise combinations of the elements of the underlying
 * range (in no defined order). A underlying range with `n` elements will therefore have `n choose 2` possible
 * combinations.
 *
 * ### Example
 *
 * \include test/snippet/range/views/pairwise_combine.cpp
 */
template <std::ranges::view underlying_range_type>
//!\cond
    requires  std::ranges::forward_range<underlying_range_type> && std::ranges::common_range<underlying_range_type>
// !\endcond
class pairwise_combine_view : public std::ranges::view_interface<pairwise_combine_view<underlying_range_type>>
{
private:

    /*!\brief The internal iterator type.
     * \tparam range_type The type of the range this iterator is operating on.
     *
     * \details
     *
     * This iterator models the iterator category of the `underlying_iterator_type`. It maintains a pair of iterators on
     * the underlying range that move over all pairwise combinations of elements. The end is reached, when the second
     * iterator points to the end of the underlying range and the first iterator to the last element of the range.
     * Also note that this iterator does not model [Cpp17Iterator](https://en.cppreference.com/w/cpp/named_req/Iterator),
     * since it does not return a reference to the represented type but a prvalue.
     * Thus this iterator might not be usable with some legacy algorithms of the STL. But it is guaranteed to work with
     * the ranges algorithms.
     */
    template <typename range_type>
    class iterator_type
    {
    private:

        //!\brief Friend declaration for iterator with different range const-ness.
        template <typename other_range_type>
            requires std::same_as<std::remove_const_t<range_type>, std::remove_const_t<other_range_type>>
        friend class iterator_type;

        //!\brief Alias type for the iterator over the passed range type.
        using underlying_iterator_type = std::ranges::iterator_t<range_type>;
        //!\brief Alias for the value type of the underlying iterator type.
        using underlying_val_t = typename std::iterator_traits<underlying_iterator_type>::value_type;
        //!\brief Alias for the reference type of the underlying iterator type.
        using underlying_ref_t = typename std::iterator_traits<underlying_iterator_type>::reference;

    public:
        /*!\name Associated types
         * \{
         */

        //!\brief The difference type.
        using difference_type   = std::ptrdiff_t;
        //!\brief The value type.
        using value_type        = std::tuple<underlying_val_t, underlying_val_t>;
        //!\brief The reference type.
        using reference         = common_tuple<underlying_ref_t, underlying_ref_t>;
        //!\brief The pointer type.
        using pointer           = void;
        //!\brief The iterator category tag.
        using iterator_category = iterator_tag_t<underlying_iterator_type>;
        //!\}

        /*!\name Constructors, destructor and assignment
         * \{
         */
        constexpr iterator_type()                                  noexcept = default; //!< Defaulted.
        constexpr iterator_type(iterator_type const &)             noexcept = default; //!< Defaulted.
        constexpr iterator_type(iterator_type &&)                  noexcept = default; //!< Defaulted.
        constexpr iterator_type & operator=(iterator_type const &) noexcept = default; //!< Defaulted.
        constexpr iterator_type & operator=(iterator_type &&)      noexcept = default; //!< Defaulted.
        ~iterator_type()                                           noexcept = default; //!< Defaulted.

        /*!\brief Constructs the iterator from the current underlying iterator and the end iterator of the underlying
         *        range.
         * \param[in] iter     The iterator pointing to current element within the underlying range.
         * \param[in] begin_it The iterator pointing to begin of the underlying range. Only needed if
         *                     `underlying_iterator_type` models std::random_access_iterator.
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

        /*!\brief Constructs const iterator from non-const iterator.
         * \param[in] other The non-const iterator to construct from.
         *
         * \details
         *
         * Allows construction of a const iterator (operating on a const range) from a non-const iterator.
         * This special constructor is needed as it is not covered by the standard copy and move constructors.
         */
        template <typename other_range_type>
        //!\cond
            requires std::convertible_to<other_range_type, range_type &> &&
                     std::same_as<std::remove_const_t<other_range_type>, std::remove_const_t<range_type>>
        //!\endcond
        constexpr iterator_type(iterator_type<other_range_type> other) noexcept :
            iterator_type{std::move(other.first_it), std::move(other.begin_it), std::move(other.end_it)}
        {}
        //!\}

        /*!\name Accessors
         * \{
         */
         //!\brief Accesses the pointed-to element
        constexpr reference operator*() const
            noexcept(noexcept(*std::declval<underlying_iterator_type>()))
        {
            return reference{*first_it, *second_it};
        }

        /*!\brief Access the element at the given index
         * \param[in] index The index of the element to be returned.
         */
        constexpr reference operator[](size_t const index)
            noexcept(noexcept(std::declval<iterator_type &>().from_index(1)))
        //!\cond
            requires std::random_access_iterator<underlying_iterator_type>
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

        //!\brief Pre-decrement operator; `underlying_iterator_type` must model std::bidirectional_iterator.
        constexpr iterator_type & operator--(/*pre-decrement*/)
            noexcept(noexcept(--std::declval<underlying_iterator_type &>()))
        //!\cond
            requires std::bidirectional_iterator<underlying_iterator_type>
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

        //!\brief Post-decrement operator; `underlying_iterator_type` must model std::bidirectional_iterator.
        constexpr iterator_type operator--(int /*post-decrement*/)
            noexcept(noexcept(std::declval<underlying_iterator_type &>()--))
        //!\cond
            requires std::bidirectional_iterator<underlying_iterator_type>
        //!\endcond
        {
            iterator_type tmp{*this};
            --*this;
            return tmp;
        }

        //!\brief Advances the iterator by the given offset; `underlying_iterator_type` must model
        //!\      std::random_access_iterator.
        constexpr iterator_type & operator+=(difference_type const offset)
            noexcept(noexcept(std::declval<iterator_type &>().from_index(1)))
        //!\cond
            requires std::random_access_iterator<underlying_iterator_type>
        //!\endcond
        {
            from_index(to_index() + offset);
            return *this;
        }

        //!\brief Advances the iterator by the given offset; `underlying_iterator_type` must model
        //!\      std::random_access_iterator.
        constexpr iterator_type operator+(difference_type const offset)
            noexcept(noexcept(std::declval<iterator_type &>() += 1))
        //!\cond
            requires std::random_access_iterator<underlying_iterator_type>
        //!\endcond
        {
            iterator_type tmp{*this};
            return (tmp += offset);
        }

        //!\brief Advances the iterator by the given offset; `underlying_iterator_type` must model
        //!\      std::random_access_iterator.
        constexpr friend iterator_type operator+(difference_type const offset, iterator_type iter)
            noexcept(noexcept(std::declval<iterator_type<range_type> &>().from_index(1)))
        //!\cond
            requires std::random_access_iterator<underlying_iterator_type>
        //!\endcond
        {
            iter.from_index(iter.to_index() + offset);
            return iter;
        }

        //!\brief Decrements the iterator by the given offset; `underlying_iterator_type` must model
        //!\      std::random_access_iterator.
        constexpr iterator_type & operator-=(difference_type const offset)
            noexcept(noexcept(std::declval<iterator_type &>().from_index(1)))
        //!\cond
            requires std::random_access_iterator<underlying_iterator_type>
        //!\endcond
        {
            from_index(to_index() - offset);
            return *this;
        }

        //!\brief Decrements the iterator by the given offset; `underlying_iterator_type` must model
        //!\      std::random_access_iterator.
        constexpr iterator_type operator-(difference_type const offset)
            noexcept(noexcept(std::declval<iterator_type &>() -= 1))
        //!\cond
            requires std::random_access_iterator<underlying_iterator_type>
        //!\endcond
        {
            iterator_type tmp{*this};
            return (tmp -= offset);
        }

        //!\brief Computes the distance between two iterators; `underlying_iterator_type` must model
        //!\      std::random_access_iterator.
        template <typename other_range_type>
        //!\cond
            requires std::random_access_iterator<underlying_iterator_type> &&
                     std::same_as<std::remove_const_t<range_type>, std::remove_const_t<other_range_type>>
        //!\endcond
        constexpr difference_type operator-(iterator_type<other_range_type> const & rhs) const
            noexcept(noexcept(std::declval<iterator_type &>().to_index()))
        {
            return static_cast<difference_type>(to_index() - rhs.to_index());
        }
        //!\}

        /*!\name Comparison operators
         * \brief These operators are available if the `underlying_iterator_type` models std::equality_comparable or
         *        std::totally_ordered respectively.
         * \{
         */
        //NOTE: The comparison operators should be implemented as friends, but due to a bug in gcc friend function
        // cannot yet be constrained. To avoid unexpected errors with the comparison all operators are implemented as
        // direct members and not as friends.

        //!\brief Checks whether `*this` is equal to `rhs`.
        template <typename other_range_type>
        //!\cond
            requires std::equality_comparable_with<underlying_iterator_type, std::ranges::iterator_t<other_range_type>> &&
                     std::same_as<std::remove_const_t<range_type>, std::remove_const_t<other_range_type>>
        //!\endcond
        constexpr bool operator==(iterator_type<other_range_type> const & rhs) const
            noexcept(noexcept(std::declval<underlying_iterator_type &>() == std::declval<underlying_iterator_type &>()))
        {
            return std::tie(first_it, second_it) == std::tie(rhs.first_it, rhs.second_it);
        }

        //!\brief Checks whether `*this` is not equal to `rhs`.
        template <typename other_range_type>
        //!\cond
            requires std::equality_comparable_with<underlying_iterator_type, std::ranges::iterator_t<other_range_type>> &&
                     std::same_as<std::remove_const_t<range_type>, std::remove_const_t<other_range_type>>
        //!\endcond
        constexpr bool operator!=(iterator_type<other_range_type> const & rhs) const
            noexcept(noexcept(std::declval<underlying_iterator_type &>() != std::declval<underlying_iterator_type &>()))
        {
            return !(*this == rhs);
        }

        //!\brief Checks whether `*this` is less than `rhs`.
        template <typename other_range_type>
        //!\cond
            requires std::totally_ordered_with<underlying_iterator_type,
                                                   std::ranges::iterator_t<other_range_type>> &&
                     std::same_as<std::remove_const_t<range_type>, std::remove_const_t<other_range_type>>
        //!\endcond
        constexpr bool operator<(iterator_type<other_range_type> const & rhs) const
            noexcept(noexcept(std::declval<underlying_iterator_type &>() < std::declval<underlying_iterator_type &>()))
        {
            return std::tie(first_it, second_it) < std::tie(rhs.first_it, rhs.second_it);
        }

        //!\brief Checks whether `*this` is greater than `rhs`.
        template <typename other_range_type>
        //!\cond
            requires std::totally_ordered_with<underlying_iterator_type,
                                                   std::ranges::iterator_t<other_range_type>> &&
                     std::same_as<std::remove_const_t<range_type>, std::remove_const_t<other_range_type>>
        //!\endcond
        constexpr bool operator>(iterator_type<other_range_type> const & rhs) const
            noexcept(noexcept(std::declval<underlying_iterator_type &>() > std::declval<underlying_iterator_type &>()))

        {
            return std::tie(first_it, second_it) > std::tie(rhs.first_it, rhs.second_it);
        }

        //!\brief Checks whether `*this` is less than or equal to `rhs`.
        template <typename other_range_type>
        //!\cond
            requires std::totally_ordered_with<underlying_iterator_type,
                                                   std::ranges::iterator_t<other_range_type>> &&
                     std::same_as<std::remove_const_t<range_type>, std::remove_const_t<other_range_type>>
        //!\endcond
        constexpr bool operator<=(iterator_type<other_range_type> const & rhs) const
            noexcept(noexcept(std::declval<underlying_iterator_type &>() <= std::declval<underlying_iterator_type &>()))
        {
            return std::tie(first_it, second_it) <= std::tie(rhs.first_it, rhs.second_it);
        }

        //!\brief Checks whether `*this` is greater than or equal to `rhs`.
        template <typename other_range_type>
        //!\cond
            requires std::totally_ordered_with<underlying_iterator_type,
                                                   std::ranges::iterator_t<other_range_type>> &&
                     std::same_as<std::remove_const_t<range_type>, std::remove_const_t<other_range_type>>
        //!\endcond
        constexpr bool operator>=(iterator_type<other_range_type> const & rhs) const
            noexcept(noexcept(std::declval<underlying_iterator_type &>() >= std::declval<underlying_iterator_type &>()))
        {
            return std::tie(first_it, second_it) >= std::tie(rhs.first_it, rhs.second_it);
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
            requires std::random_access_iterator<underlying_iterator_type>
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
            requires std::random_access_iterator<underlying_iterator_type>
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

    /*!\name Associated types
     * \{
     */

    //!\brief The iterator type.
    using iterator          = iterator_type<underlying_range_type>;
    //!\brief The const iterator type. Evaluates to void if the underlying range is not const iterable.
    using const_iterator    = transformation_trait_or_t<std::type_identity<iterator_type<underlying_range_type const>>,
                                                        void>;
    //!\brief The reference_type.
    using reference         = typename iterator::reference;
    //!\brief The const_reference type is equal to the reference type.
    using const_reference   = reference;
    //!\brief The value_type (which equals the reference_type with any references removed).
    using value_type        = typename iterator::value_type;
    //!\brief If the underliying range is Sized, this resolves to range_type::size_type, otherwise void.
    using size_type         = detail::transformation_trait_or_t<seqan3::size_type<underlying_range_type>, void>;
    //!\brief A signed integer type, usually std::ptrdiff_t.
    using difference_type   = std::iter_difference_t<iterator>;
    //!\}

    /*!\name Constructors, destructor and assignment
     * \{
     */
    //!\brief Default Default-Constructor.
    constexpr pairwise_combine_view() = default;
    //!\brief Default Copy-Constructor.
    constexpr pairwise_combine_view(pairwise_combine_view const &) = default;
    //!\brief Default Move-Constructor.
    constexpr pairwise_combine_view(pairwise_combine_view &&) = default;
    //!\brief Default Copy-Assignment.
    constexpr pairwise_combine_view & operator=(pairwise_combine_view const &) = default;
    //!\brief Default Move-Assignment.
    constexpr pairwise_combine_view & operator=(pairwise_combine_view &&) = default;
    //!\brief Default Destructor.
    ~pairwise_combine_view() = default;

    /*!\brief Constructs from a view.
     * \param[in] range The underlying range to be wrapped. Of type `underlying_range_type`.
     *
     * \details
     *
     * During construction the iterator pointing to the last element of the view is cached (not the end of the range).
     * This optimises the call to end if the underlying range models only forward_range. Otherwise the call to end will be
     * linear in the number of elements of the underlying range.
     *
     * \attention This view cannot be chained immediately after an infinite range, because upon construction it will
     *            take forever to reach the last element of the view.
     *
     * ### Complexity
     *
     * Constant if `underlying_range_type` models std::ranges::bidirectional_range, otherwise linear.
     */
    explicit constexpr pairwise_combine_view(underlying_range_type range) : u_range{std::move(range)}
    {
        // Check if range is empty.
        if (std::ranges::empty(u_range))
        {
            back_iterator = std::ranges::end(u_range);
        }
        else
        {
            if constexpr (std::ranges::bidirectional_range<underlying_range_type>)
            { // Simply take one before the end. We can do this as we require underlying_range_type to be a common range.
                back_iterator = std::ranges::prev(std::ranges::end(u_range));
            }
            else
            { // For all other cases we need to set the back_iterator in linear time to the correct position.
                back_iterator = std::ranges::begin(u_range);
                if constexpr (std::ranges::sized_range<underlying_range_type>)
                {
                    std::ranges::advance(back_iterator, std::ranges::size(u_range) - 1);
                }
                else // We don't have the size, so we need to increment until one before the end in a linear pass.
                {
                    auto tmp_it = back_iterator;
                    do
                    {
                        back_iterator = tmp_it;
                    } while (++tmp_it != std::ranges::end(u_range));
                }
            }
        }
    }

    /*!\brief Constructs from a view.
     * \tparam    other_range_t  The type of the range to be wrapped with seqan3::detail::pairwise_combine_view;
     *                           must model std::ranges::viewable_range and underlying_range_type must be constructible
     *                           with other_range wrapped in std::views::all.
     * \param[in] range          The underlying range to be wrapped.
     *
     * \details
     *
     * During construction the iterator pointing to the last element of the view is cached (not the end of the range).
     * This optimises the call to end if the underlying range models only forward_range. Otherwise the call to end will be
     * linear in the number of elements of the underlying range.
     *
     * \attention This view cannot be chained immediately after an infinite range, because upon construction it will
     *            take forever to reach the last element of the view.
     *
     * ### Complexity
     *
     * Constant if `other_range_t` models std::ranges::bidirectional_range, otherwise linear.
     */
    template <typename other_range_t>
    //!\cond
        requires !std::same_as<remove_cvref_t<other_range_t>, pairwise_combine_view> &&
                 std::ranges::viewable_range<other_range_t> &&  // Must come after self type check to avoid conflicts with the move constructor.
                 std::constructible_from<underlying_range_type, ranges::ref_view<std::remove_reference_t<other_range_t>>>
            //TODO: Investigate: the following expression is equivalent to the one above but raises a weird assertion in
            //      the ranges adaptor suggesting that the pairwise_combine_view is not a viewable_range.
            //      std::constructible_from<underlying_range_type, decltype(std::views::all(std::declval<other_range_t &&>()))>
    //!\endcond
    explicit constexpr pairwise_combine_view(other_range_t && range) :
        pairwise_combine_view{std::views::all(std::forward<other_range_t>(range))}
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
    constexpr iterator begin() noexcept
    {
        return {std::ranges::begin(u_range), std::ranges::begin(u_range), std::ranges::end(u_range)};
    }

    //!\copydoc begin()
    constexpr const_iterator begin() const noexcept
    //!\cond
        requires const_iterable_range<underlying_range_type>
    //!\endcond
    {
        return {std::ranges::begin(u_range), std::ranges::begin(u_range), std::ranges::end(u_range)};
    }

    //!\copydoc begin()
    constexpr const_iterator cbegin() const noexcept
    //!\cond
        requires const_iterable_range<underlying_range_type>
    //!\endcond
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
    constexpr iterator end() noexcept
    {
        return {back_iterator, std::ranges::begin(u_range), std::ranges::end(u_range)};
    }

    //!\copydoc end()
    constexpr const_iterator end() const noexcept
    //!\cond
        requires const_iterable_range<underlying_range_type>
    //!\endcond
    {
        return {back_iterator, std::ranges::begin(u_range), std::ranges::end(u_range)};
    }

    //!\copydoc end()
    constexpr const_iterator cend() const noexcept
    //!\cond
        requires const_iterable_range<underlying_range_type>
    //!\endcond
    {
        return end();
    }
    //!\}

    /*!\name Capacity
     * \{
     */
    //!\brief Computes the size based on the size of the underlying range.
    constexpr size_type size() const noexcept
    //!\cond
        requires std::ranges::sized_range<underlying_range_type>
    //!\endcond
    {
        return (std::ranges::size(u_range) * (std::ranges::size(u_range) - 1) / 2);
    }
    //!\}

private:

    //!\brief The underling range.
    underlying_range_type u_range{};
    //!\brief The cached iterator pointing to the last element of the underlying range.
    std::ranges::iterator_t<underlying_range_type> back_iterator{};
};

/*!\name Type deduction guides
 * \{
 */

//!\brief Deduces the correct template type from a non-view lvalue range by wrapping the range in std::views::all.
template <std::ranges::viewable_range other_range_t>
pairwise_combine_view(other_range_t && range) ->
    pairwise_combine_view<std::ranges::all_view<other_range_t>>;
//!\}

} // namespace seqan3::detail

namespace seqan3::views
{
/*!\name General purpose views
 * \{
 */

/*!\brief             A view adaptor that generates all pairwise combinations of the elements of the underlying range.
 * \tparam urng_t     The type of the range being processed. See below for requirements.
 * \param[in] urange  The range being processed.
 * \returns           A view over all pairwise combinations of the elements of the underlying range.
 *                    See below for the properties of the returned range.
 * \ingroup views
 *
 * \details
 *
 * \header_file{seqan3/range/views/pairwise_combine.hpp}
 *
 * This view generates two-element tuples representing all unique combinations of the elements of the underlying range
 * (the order of the elements is undefined). If the underlying range has less than two elements the returned range is
 * empty, otherwise the size of the returned range corresponds to the binomial coefficient `n choose 2`, where `n` is
 * the size of the underlying range. The reference type of this range is a tuple over the reference type of the
 * underlying range.
 * In order to receive the end iterator in constant time an iterator pointing to the last element of the underlying
 * range will be cached upon construction of this view. This construction takes linear time for underlying ranges that
 * do not model std::ranges::bidirectional_range.
 *
 * ### Iterator
 *
 * The returned iterator from begin does not model [Cpp17Iterator](https://en.cppreference.com/w/cpp/named_req/Iterator)
 * since it does not return a reference to the represented type but a prvalue. Thus this iterator might not be usable
 * within some legacy algorithms of the STL. But it is guaranteed to work with the ranges algorithms.
 *
 * ### View properties
 *
 * | Concepts and traits              | `urng_t` (underlying range type)      | `rrng_t` (returned range type)                                                               |
 * |----------------------------------|:-------------------------------------:|:--------------------------------------------------------------------------------------------:|
 * | std::ranges::input_range         | *required*                            | *preserved*                                                                                  |
 * | std::ranges::forward_range       | *required*                            | *preserved*                                                                                  |
 * | std::ranges::bidirectional_range |                                       | *preserved*                                                                                  |
 * | std::ranges::random_access_range |                                       | *preserved*                                                                                  |
 * | std::ranges::contiguous_range    |                                       | *lost*                                                                                       |
 * |                                  |                                       |                                                                                              |
 * | std::ranges::viewable_range      | *required*                            | *guaranteed*                                                                                 |
 * | std::ranges::view                |                                       | *guaranteed*                                                                                 |
 * | std::ranges::sized_range         |                                       | *preserved*                                                                                  |
 * | std::ranges::common_range        | *required*                            | *guaranteed*                                                                                 |
 * | std::ranges::output_range        |                                       | *preserved*                                                                                  |
 * | seqan3::const_iterable_range     |                                       | *preserved*                                                                                  |
 * |                                  |                                       |                                                                                              |
 * | std::ranges::range_reference_t   |                                       | common_tuple<std::ranges::range_reference_t<urng_t>, std::ranges::range_reference_t<urng_t>> |
 *
 * See the \link views views submodule documentation \endlink for detailed descriptions of the view properties.
 *
 * ### Thread safety
 *
 * Concurrent access to this view, e.g. while iterating over it, is thread-safe and must not be protected externally.
 *
 * ### Example
 *
 * \include test/snippet/range/views/pairwise_combine.cpp
 *
 * \attention This view cannot be chained immediately with an infinite range, because upon construction it will
 *            take forever to reach the last element of the view.
 *
 * \hideinitializer
 */
inline constexpr auto pairwise_combine = detail::adaptor_for_view_without_args<detail::pairwise_combine_view>{};

//!\}
} // namespace seqan3::views
