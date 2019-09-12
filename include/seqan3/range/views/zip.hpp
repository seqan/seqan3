// -----------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2020, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2020, Knut Reinert & MPI für molekulare Genetik
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at: https://github.com/seqan/seqan3/blob/master/LICENSE.md
// -----------------------------------------------------------------------------------------------------

/*!\file
 * \brief Provides seqan3::views::zip.
 * \author Enrico Seiler <enrico.seiler AT fu-berlin.de>
 */

#pragma once

#include <seqan3/core/common_tuple.hpp>
#include <seqan3/range/concept.hpp>
#include <seqan3/range/views/detail.hpp>
#include <seqan3/std/ranges>

//!\cond
namespace seqan3::detail
{

template <typename t>
SEQAN3_CONCEPT simple_view = std::ranges::view<t> &&
                             std::ranges::range<t const> &&
                             std::same_as<std::ranges::iterator_t<t>, std::ranges::iterator_t<t const>> &&
                             std::same_as<std::ranges::sentinel_t<t>, std::ranges::sentinel_t<t const>>;

template <std::weakly_incrementable ...t>
struct common_difference
{
    using type = std::conditional_t<sizeof...(t) == 0,
                                    int,
                                    std::common_type_t<std::ranges::iter_difference_t<t>...>>;
};

template <std::weakly_incrementable ...t>
using common_difference_t = typename common_difference<t...>::type;
//!\endcond

// ---------------------------------------------------------------------------------------------------------------------
// zip_view class
// ---------------------------------------------------------------------------------------------------------------------

/*!\brief The type returned by seqan3::views::zip.
 * \tparam urng_t The types of the underlying ranges, must satisfy std::ranges::view.
 * \implements std::ranges::view
 * \implements std::ranges::random_access_range
 * \implements std::ranges::sized_range
 * \implements std::ranges::output_range
 * \ingroup views
 *
 * \details
 *
 * Note that most members of this class are generated by ranges::view_interface which is not yet documented here.
 */
template <std::ranges::view ...urng_t>
class zip_view : public std::ranges::view_interface<zip_view<urng_t...>>
{
private:
    //!\cond
    template <bool constness, typename indices>
    struct zip_iterator;

    template <bool constness, typename indices>
    struct zip_sentinel;
    //!\endcond

    //!\brief The type of the std::index_sequence for `urng_t`.
    using indices = std::index_sequence_for<urng_t...>;

    //!\brief Indicates whether all types in `urng_t` model simple_view.
    static constexpr bool all_simple = (simple_view<urng_t> && ...);
    //!\brief Indicates whether all types in `urng_t` model std::ranges::forward_range.
    static constexpr bool all_forward = (std::ranges::forward_range<urng_t> && ...);
    //!\brief Indicates whether all types in `urng_t const` model std::ranges::forward_range.
    static constexpr bool all_forward_const = (std::ranges::forward_range<urng_t const> && ...);
    //!\brief Indicates whether all types in `urng_t` model std::ranges::sized_range.
    static constexpr bool all_sized = (std::ranges::sized_range<urng_t> && ...);
    //!\brief Indicates whether all types in `urng_t const` model std::ranges::sized_range.
    static constexpr bool all_sized_const = (std::ranges::sized_range<urng_t const> && ...);

    //!\brief The underlying ranges.
    std::tuple<urng_t...> urange{};
    //!\brief The iterators of the underlying ranges.
    std::tuple<std::ranges::iterator_t<urng_t>...> iterators{};

    /*!\brief Helper function for calling begin on each iterator in `iterators`.
     * \tparam N Pack of size_t.
     * \param index_seq The std::index_sequence over `N`.
     */
    template <size_t ...N>
    constexpr void begin_impl(std::index_sequence<N...> SEQAN3_DOXYGEN_ONLY(index_seq)) noexcept
    {
        ((std::get<N>(iterators) = std::ranges::begin(std::get<N>(urange))), ...);
    }

    /*!\brief Helper function for determining the smallest range in `urng_t`.
     * \tparam N Pack of size_t.
     * \param index_seq The std::index_sequence over `N`.
     */
    template <size_t ...N>
    constexpr auto size_impl(std::index_sequence<N...> SEQAN3_DOXYGEN_ONLY(index_seq)) const noexcept
    {
        if constexpr(sizeof...(urng_t) == 0)
            return 0;
        else
            return static_cast<common_difference_t<std::ranges::iterator_t<urng_t>...>>(
                       std::min({std::ranges::size(std::get<N>(urange))...}));
    }

public:
    /*!\name Associated types
     * \{
     */
    //!\brief The iterator type.
    using iterator          = zip_iterator<false, indices>;
    //!\brief The const iterator type.
    using const_iterator    = zip_iterator<true, indices>;
    //!\brief The reference_type.
    using reference         = typename iterator::reference;
    //!\brief The const_reference type.
    using const_reference   = typename const_iterator::reference;
    //!\brief The value_type.
    using value_type        = typename iterator::value_type;
    //!\brief The size type.
    using size_type         = size_t;
    //!\brief A signed integer type, usually std::ptrdiff_t.
    using difference_type   = typename iterator::difference_type;
    //!\}

    /*!\name Constructors, destructor and assignment
     * \{
     */
    constexpr zip_view()                                 = default; //!< Defaulted.
    constexpr zip_view(zip_view const & rhs)             = default; //!< Defaulted.
    constexpr zip_view(zip_view && rhs)                  = default; //!< Defaulted.
    constexpr zip_view & operator=(zip_view const & rhs) = default; //!< Defaulted.
    constexpr zip_view & operator=(zip_view && rhs)      = default; //!< Defaulted.
    ~zip_view()                                          = default; //!< Defaulted.

    //!\brief Construct from another view.
    constexpr zip_view(urng_t... rng) noexcept((std::is_nothrow_move_constructible_v<urng_t> && ...))
    //!\cond
        requires sizeof...(urng_t) != 0
    //!\endcond
    :urange{std::move(rng)...}
    {}
    //!\}

    /*!\name Iterators
     * \{
     */
    /*!\brief Returns an iterator to the first element of the container.
     * \returns Iterator to the first element.
     *
     * If the ranges are empty, the returned iterator will be equal to end().
     *
     * ### Complexity
     *
     * Constant.
     *
     * ### Exceptions
     *
     * No-throw guarantee.
     */
    constexpr auto begin() noexcept
    {
        if constexpr (!all_forward)
            begin_impl(indices{});

        return zip_iterator<all_simple && all_forward, indices>{*this};
    }

    //!\copydoc begin()
    constexpr auto begin() const noexcept
    //!\cond
        requires all_forward_const && (const_iterable_range<urng_t> && ...)
    //!\endcond
    {
        return zip_iterator<true, indices>{*this};
    }

    //!\copydoc begin()
    constexpr auto cbegin() const noexcept
    //!\cond
        requires all_forward_const && (const_iterable_range<urng_t> && ...)
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
    constexpr auto end() noexcept
    //!\cond
        requires all_forward
    //!\endcond
    {
        if constexpr((std::ranges::random_access_range<urng_t> && ...) && all_sized)
            return zip_iterator<all_simple, indices>{*this, size_impl(indices{})};
        else
            return zip_sentinel<all_simple, indices>{*this};
    }

    //!\copydoc end()
    constexpr auto end() const noexcept
    //!\cond
        requires all_forward && (const_iterable_range<urng_t> && ...)
    //!\endcond
    {
        if constexpr((std::ranges::random_access_range<urng_t const> && ...) && all_sized_const)
            return zip_iterator<true, indices>{*this, size_impl(indices{})};
        else if constexpr(all_forward && all_forward_const)
            return zip_sentinel<all_simple, indices>{*this};
        else
            return std::ranges::default_sentinel;
    }

    //!\copydoc end()
    constexpr auto cend() const noexcept
    //!\cond
        requires all_forward && (const_iterable_range<urng_t> && ...)
    //!\endcond
    {
        return end();
    }
    //!\}

    /*!\brief Returns the number of elements of the smallest range.
     * \returns The number of elements of the smallest range.
     *
     * ### Complexity
     *
     * Constant.
     *
     * ### Exceptions
     *
     * No-throw guarantee.
     */
    constexpr auto size() noexcept
    //!\cond
        requires all_sized
    //!\endcond
    {
        return size_impl(indices{});
    }

    //!\copydoc size()
    constexpr auto size() const noexcept
    //!\cond
        requires all_sized_const
    //!\endcond
    {
        return size_impl(indices{});
    }
};

//!\brief Template argument type deduction guide that wraps ranges in std::views::all.
template <std::ranges::viewable_range ...rng_t>
zip_view(rng_t&&...) -> zip_view<std::ranges::all_view<rng_t>...>;

/*!\brief The iterator type of seqan3::views::zip.
 * \tparam constness Indicates whether the iterator is const.
 * \tparam N Pack of size_t.
 */
template <std::ranges::view ...urng_t>
template <bool constness, size_t ...N>
class zip_view<urng_t...>::zip_iterator<constness, std::index_sequence<N...>>
{
private:
    //!\brief Helper function that returns `t` or `t const` depending on `constness`.
    template <typename t>
    using maybe_const = std::conditional_t<constness, t const, t>;

    //!\brief The type of the parent view.
    using parent_t = maybe_const<zip_view<urng_t...>>;

    //!\brief The type of the sentinel.
    using sentinel = typename zip_view<urng_t...>::template
                     zip_sentinel<constness, std::index_sequence<N...>>;

    //!\brief Indicates whether all types in `maybe_const<urng_t>` model std::ranges::random_access_range.
    static constexpr bool all_random_access = (std::ranges::random_access_range<maybe_const<urng_t>> && ...);
    //!\brief Indicates whether all types in `maybe_const<urng_t>` model std::ranges::bidirectional_range.
    static constexpr bool all_bidi = (std::ranges::bidirectional_range<maybe_const<urng_t>> && ...);
    //!\brief Indicates whether all types in `maybe_const<urng_t>` model std::ranges::forward_range.
    static constexpr bool all_forward = (std::ranges::forward_range<maybe_const<urng_t>> && ...);
    //!\brief Indicates whether all types in `maybe_const<urng_t>` model std::ranges::input_range.
    static constexpr bool all_input = (std::ranges::input_range<maybe_const<urng_t>> && ...);

    static_assert(!constness || all_forward,
                  "The zip_iterator must either be not const or all ranges must model std::ranges::forward_range.");
    static_assert(std::same_as<std::index_sequence_for<urng_t...>, std::index_sequence<N...>>,
                  "The number of ranges differes from the passed template parameter.");

    //!\brief A pointer to the parent view.
    parent_t * parent = nullptr;
    //!\brief The iterators of each range in `urng_t`.
    std::tuple<std::ranges::iterator_t<maybe_const<urng_t>>...> iterators{};

public:
    /*!\name Associated types
     * \{
     */
    //!\brief Type for distances between iterators.
    using difference_type = common_difference_t<std::ranges::iterator_t<maybe_const<urng_t>>...>;
    //!\brief Value type of this iterator.
    using value_type = std::tuple<std::ranges::iter_value_t<std::ranges::iterator_t<maybe_const<urng_t>>>...>;
    //!\brief The pointer type.
    using pointer = void;
    //!\brief Reference to `value_type`.
    using reference = common_tuple<std::ranges::iter_reference_t<std::ranges::iterator_t<maybe_const<urng_t>>>...>;
    //!\brief Reference to `value_type`.
    using const_reference = common_tuple<std::ranges::iter_reference_t<std::ranges::iterator_t<
                                maybe_const<urng_t>>> const...>;
    //!\brief Tag this class depending on which concept are modelled by `urng_t`.
    using iterator_category = std::conditional_t<all_random_access, std::random_access_iterator_tag,
                                std::conditional_t<all_bidi, std::bidirectional_iterator_tag,
                                    std::conditional_t<all_forward, std::forward_iterator_tag,
                                        std::conditional_t<all_input, std::input_iterator_tag,
                                            std::output_iterator_tag>>>>;
    //!\}

    /*!\name Constructors, destructor and assignment
     * \{
     */
    constexpr zip_iterator()                                 = default; //!< Defaulted.
    constexpr zip_iterator(zip_iterator const &)             = default; //!< Defaulted.
    constexpr zip_iterator(zip_iterator &&)                  = default; //!< Defaulted.
    constexpr zip_iterator & operator=(zip_iterator const &) = default; //!< Defaulted.
    constexpr zip_iterator & operator=(zip_iterator &&)      = default; //!< Defaulted.
    ~zip_iterator()                                          = default; //!< Defaulted.

    //!\brief Construct from reference to view.
    explicit constexpr zip_iterator(parent_t & parent_) : parent{std::addressof(parent_)} {};

    //!\brief Construct from reference to view.
    //!       Sets all iterators to respective `begin` if all ranges model std::ranges::forward_range.
    constexpr zip_iterator(parent_t & parent_)
    //!\cond
        requires all_forward
    //!\endcond
    : parent{std::addressof(parent_)}, iterators{std::ranges::begin(std::get<N>(parent->urange)) ...} {};

    //!\brief Construct from reference to view and an offset.
    //!       Sets all iterators to the nth position if all ranges model std::ranges::random_access_range.
    constexpr zip_iterator(parent_t & parent_, difference_type n)
    //!\cond
        requires all_random_access
    //!\endcond
    : parent{std::addressof(parent_)}, iterators{(std::ranges::begin(std::get<N>(parent->urange)) + n) ...} {};
    //!\}

    //!\brief Dereferences the iterators.
    constexpr auto operator*() noexcept(noexcept((*std::declval<std::ranges::iterator_t<urng_t>>(), ...)))
    {
        return reference{*std::get<N>(iterators)...};
    }

    //!\brief Dereferences the iterators.
    constexpr auto operator*() const noexcept(noexcept((*std::declval<std::ranges::iterator_t<urng_t const>>(), ...)))
    //!\cond
        requires (std::input_or_output_iterator<std::ranges::iterator_t<maybe_const<urng_t const>>> && ...)
    //!\endcond
    {
        return const_reference{*std::get<N>(iterators)...};
    }

    /*!\name Arithmetic operators
     * \{
     */

    //!\brief Increments the iterator by one.
    constexpr zip_iterator& operator++() noexcept(noexcept((++std::declval<std::ranges::iterator_t<urng_t>>(), ...)))
    {
        (++std::get<N>(iterators), ...);
        return *this;
    }

    //!\brief Returns an iterator incremented by one.
    constexpr auto operator++(int) noexcept(noexcept((++std::declval<std::ranges::iterator_t<urng_t>>(), ...)) &&
                                       (std::is_nothrow_copy_constructible_v<std::ranges::iterator_t<urng_t>> && ...))
    {
        zip_iterator tmp{*this};
        ++*this;
        return tmp;
    }

    //!\brief Decrements the iterator by one.
    constexpr zip_iterator& operator--() noexcept(noexcept((--std::declval<std::ranges::iterator_t<urng_t>>(), ...)))
    //!\cond
        requires all_bidi
    //!\endcond
    {
        (--std::get<N>(iterators), ...);
        return *this;
    }

    //!\brief Returns an iterator decremented by one.
    constexpr zip_iterator operator--(int) noexcept(noexcept(
        (--std::declval<std::ranges::iterator_t<urng_t>>(), ...)) &&
        (std::is_nothrow_copy_constructible_v<std::ranges::iterator_t<urng_t>> && ...))
    //!\cond
        requires all_bidi
    //!\endcond
    {
        zip_iterator tmp{*this};
        --*this;
        return tmp;
    }

    //!\brief Advances the iterator by `n` positions.
    constexpr zip_iterator& operator+=(difference_type n) noexcept(noexcept(
        ((std::declval<std::ranges::iterator_t<urng_t>>() += n), ...)))
    //!\cond
        requires all_random_access
    //!\endcond
    {
        ((std::get<N>(iterators) += n), ...);
        return *this;
    }

    //!\brief Returns an iterator advanced by `n` positions.
    friend constexpr zip_iterator operator+(zip_iterator const & lhs, difference_type n) noexcept(noexcept(
        ((std::declval<std::ranges::iterator_t<urng_t>>() += n), ...)) &&
        (std::is_nothrow_copy_constructible_v<std::ranges::iterator_t<urng_t>> && ...))
    //!\cond
        requires all_random_access
    //!\endcond
    {
        zip_iterator tmp{lhs};
        return tmp += n;
    }

    //!\brief Returns an iterator advanced by `n` positions.
    friend constexpr zip_iterator operator+(difference_type n, zip_iterator const & rhs) noexcept(noexcept(
        ((std::declval<std::ranges::iterator_t<urng_t>>() += n), ...)) &&
        (std::is_nothrow_copy_constructible_v<std::ranges::iterator_t<urng_t>> && ...))
    //!\cond
        requires all_random_access
    //!\endcond
    {
        return rhs + n;
    }

    //!\brief Advances the iterator by `-n` positions.
    constexpr zip_iterator& operator-=(difference_type n) noexcept(noexcept(
        ((std::declval<std::ranges::iterator_t<urng_t>>() -= n), ...)))
    //!\cond
        requires all_random_access
    //!\endcond
    {
        ((std::get<N>(iterators) += n), ...);
        return *this;
    }

    //!\brief Returns an iterator advanced by `-n` positions.
    friend constexpr zip_iterator operator-(zip_iterator const & lhs, difference_type n) noexcept(noexcept(
        ((std::declval<std::ranges::iterator_t<urng_t>>() -= n), ...)) &&
        (std::is_nothrow_copy_constructible_v<std::ranges::iterator_t<urng_t>> && ...))
    //!\cond
        requires all_random_access
    //!\endcond
    {
        zip_iterator tmp{lhs};
        return tmp -= n;
    }

    //!\brief Returns the distance between two iterators.
    friend constexpr difference_type operator-(zip_iterator const & lhs, zip_iterator const & rhs) noexcept(noexcept(
        ((std::declval<std::ranges::iterator_t<urng_t>>() - std::declval<std::ranges::iterator_t<urng_t>>()), ...)))
    //!\cond
        requires (std::sized_sentinel_for<std::ranges::iterator_t<maybe_const<urng_t>>,
                                          std::ranges::iterator_t<maybe_const<urng_t>>> && ...)
    //!\endcond
    {
        if constexpr(sizeof...(urng_t) == 0)
            return 0;
        else
            return std::max({std::abs(std::get<N>(lhs.iterators) - std::get<N>(rhs.iterators))...});
    }

    //!\brief Accesses an element by index.
    constexpr auto operator[](difference_type n) const noexcept(noexcept(
        (std::declval<std::ranges::iterator_t<urng_t>>()[n], ...)))
    //!\cond
        requires all_random_access
    //!\endcond
    {
        return *(*this + n);
    }

    //!\brief Checks whether `lhs` is equal to `rhs`.
    friend constexpr bool operator==(zip_iterator const & lhs, zip_iterator const & rhs) noexcept
    //!\cond
        requires all_forward
    //!\endcond
    {
        return lhs.iterators == rhs.iterators;
    }

    //!\brief Checks whether `lhs` is not equal to `rhs`.
    friend constexpr bool operator!=(zip_iterator const & lhs, zip_iterator const & rhs) noexcept
    //!\cond
        requires all_forward
    //!\endcond
    {
        return !(lhs == rhs);
    }

    //!\brief Checks whether any iterator of `lhs` is equal to the respective end().
    friend constexpr bool operator==(zip_iterator const & lhs, std::ranges::default_sentinel_t) noexcept
    //!\cond
        requires !all_forward
    //!\endcond
    {
        return ((std::get<N>(lhs.iterators) == std::ranges::end(std::get<N>(lhs.parent->uranges))) || ...);
    }

    //!\brief Checks whether any iterator of `rhs` is equal to the respective end().
    friend constexpr bool operator==(std::ranges::default_sentinel_t lhs, zip_iterator const & rhs) noexcept
    //!\cond
        requires !all_forward
    //!\endcond
    {
        return rhs == lhs;
    }

    //!\brief Checks whether `lhs` is unequal to `rhs`.
    friend constexpr bool operator!=(zip_iterator const & lhs, std::ranges::default_sentinel_t rhs) noexcept
    //!\cond
        requires !all_forward
    //!\endcond
    {
        return !(lhs == rhs);
    }

    //!\brief Checks whether `lhs` is unequal to `rhs`.
    friend constexpr bool operator!=(std::ranges::default_sentinel_t lhs, zip_iterator const & rhs) noexcept
    //!\cond
        requires !all_forward
    //!\endcond
    {
        return !(rhs == lhs);
    }

    //!\brief Checks whether `lhs` is equal to `rhs`.
    friend constexpr bool operator==(zip_iterator const & lhs, sentinel const & rhs) noexcept
    //!\cond
        requires all_forward
    //!\endcond
    {
        return ((std::get<N>(lhs.iterators) == std::get<N>(rhs.iterators)) || ...);
    }

    //!\brief Checks whether `lhs` is equal to `rhs`.
    friend constexpr bool operator==(sentinel const & lhs, zip_iterator const & rhs) noexcept
    //!\cond
        requires all_forward
    //!\endcond
    {
        return rhs == lhs;
    }

    //!\brief Checks whether `lhs` is unequal to `rhs`.
    friend constexpr bool operator!=(zip_iterator const & lhs, sentinel const & rhs) noexcept
    //!\cond
        requires all_forward
    //!\endcond
    {
        return !(lhs == rhs);
    }

    //!\brief Checks whether `lhs` is unequal to `rhs`.
    friend constexpr bool operator!=(sentinel const & lhs, zip_iterator const & rhs) noexcept
    //!\cond
        requires all_forward
    //!\endcond
    {
        return !(rhs == lhs);
    }

    //!\brief Checks whether `lhs` is smaller than `rhs`.
    friend constexpr bool operator<(zip_iterator const & lhs, zip_iterator const & rhs) noexcept
    //!\cond
        requires all_random_access
    //!\endcond
    {
        return lhs.iterators < rhs.iterators;
    }

    //!\brief Checks whether `lhs` is bigger than `rhs`.
    friend constexpr bool operator>(zip_iterator const & lhs, zip_iterator const & rhs) noexcept
    //!\cond
        requires all_random_access
    //!\endcond
    {
        return rhs < lhs;
    }

    //!\brief Checks whether `lhs` is smaller than or equal to `rhs`.
    friend constexpr bool operator<=(zip_iterator const & lhs, zip_iterator const & rhs) noexcept
    //!\cond
        requires all_random_access
    //!\endcond
    {
        return !(rhs < lhs);
    }

    //!\brief Checks whether `lhs` is bigger than or equal to `rhs`.
    friend constexpr bool operator>=(zip_iterator const & lhs, zip_iterator const & rhs) noexcept
    //!\cond
        requires all_random_access
    //!\endcond
    {
        return !(lhs < rhs);
    }
};

/*!\brief The sentinel type of seqan3::views::zip.
 * \tparam constness Indicates whether the iterator is const.
 * \tparam N Pack of size_t.
 */
template <std::ranges::view ...urng_t>
template <bool constness, size_t ...N>
class zip_view<urng_t...>::zip_sentinel<constness, std::index_sequence<N...>>
{
private:
    //!\brief Helper function that returns `t` or `t const` depending on `constness`.
    template <typename t>
    using maybe_const = std::conditional_t<constness, t const, t>;

    //!\brief The iterator type of the zip_view.
    using iterator = typename zip_view<urng_t...>::template
                     zip_iterator<constness, std::index_sequence<N...>>;
    //!\brief Befriend the zip_iterator.
    friend iterator;

    //!\brief The type of the parent view.
    using parent_t = zip_view<urng_t...>;

    //!\brief Indicates whether all types in `maybe_const<urng_t>` model std::ranges::forward_range.
    static constexpr bool all_forward = (std::ranges::forward_range<maybe_const<urng_t>> && ...);
    static_assert(all_forward, "The zip_sentinel requires all ranges to model std::ranges::forward_range.");
    static_assert(std::same_as<std::index_sequence_for<urng_t...>, std::index_sequence<N...>>,
                  "The number of ranges differes from the passed template parameter.");

    //!\brief The sentinels of each range in `urng_t`.
    std::tuple<std::ranges::sentinel_t<maybe_const<urng_t>>...> iterators{};

public:
    /*!\name Constructors, destructor and assignment
     * \{
     */
    constexpr zip_sentinel()                                 = default; //!< Defaulted.
    constexpr zip_sentinel(zip_sentinel const &)             = default; //!< Defaulted.
    constexpr zip_sentinel(zip_sentinel &&)                  = default; //!< Defaulted.
    constexpr zip_sentinel & operator=(zip_sentinel const &) = default; //!< Defaulted.
    constexpr zip_sentinel & operator=(zip_sentinel &&)      = default; //!< Defaulted.
    ~zip_sentinel()                                          = default; //!< Defaulted.

    //!\brief Construct from reference to view. Sets all iterators to respective `end`.
    explicit constexpr zip_sentinel(parent_t const & parent)
        : iterators{std::ranges::end(std::get<N>(parent.urange)) ...}
    {}
    //!\}
};

// ============================================================================
//  zip_fn (adaptor definition)
// ============================================================================

//!\brief View adaptor definition for views::zip.
struct zip_fn
{
    //!\brief Returns an instance of seqan3::detail::zip_view constructed with \p rngs.
    template <typename ...rng_t>
    //!\cond
        requires (std::ranges::range<rng_t> && ...)
    //!\endcond
    constexpr auto operator()(rng_t&&... rngs) const
    {
        return zip_view{std::forward<rng_t>(rngs)...};
    }
};
}

namespace seqan3::views
{

/*!\brief A range adaptor that transforms a tuple of range into a range of tuples.
 * \tparam          urng_ts The types of the ranges being processed. See below for requirements.
 * \param[in]       uranges The ranges being processed.
 * \returns A view of n-sized-tuples where n is the number of underlying ranges the i-thof size n of the elements
 *          produced by applied the invocable to each element in the underlying range.
 * \ingroup views
 *
 * \details
 *
 * **Header**
 * ```cpp
 *      #include <seqan3/range/views/zip.hpp>
 * ```
 *
 * ### View properties
 *
 * This view is **source-only**, it can only be at the beginning of a pipe of range transformations. The
 * "underlying ranges" refer to the mandatory parameters of this adaptor, not a previous range in a pipe.
 *
 * | Concepts and traits              | `urng_ts` (underlying ranges) | `rrng_t` (returned range type)             |
 * |----------------------------------|:-----------------------------:|:-------------------------------------------|
 * | std::ranges::input_range         | *required*                    | *preserved*                                |
 * | std::ranges::forward_range       |                               | *preserved*                                |
 * | std::ranges::bidirectional_range |                               | *preserved*                                |
 * | std::ranges::random_access_range |                               | *preserved*                                |
 * | std::ranges::contiguous_range    |                               | *lost*                                     |
 * |                                  |                               |                                            |
 * | std::ranges::viewable_range      | *required*                    | *guaranteed*                               |
 * | std::ranges::view                |                               | *guaranteed*                               |
 * | std::ranges::sized_range         |                               | *preserved*                                |
 * | std::ranges::common_range        |                               | *preserved*                                |
 * | std::ranges::output_range        |                               | *preserved*                                |
 * | seqan3::const_iterable_range     |                               | *preserved*                                |
 * |                                  |                               |                                            |
 * | std::ranges::range_reference_t   |                               | common_tuple<std::reference_t<urng_ts...>> |
 *
 * The guarantees for the returned range type only hold if the respective requirements are met by **all underlying
 * ranges**.
 *
 * See the \link views views submodule documentation \endlink for detailed descriptions of the view properties.
 *
 * ### Example
 *
 * \include test/snippet/range/views/zip.cpp
 *
 * \hideinitializer
 */
inline constexpr auto zip = detail::zip_fn{};

}
