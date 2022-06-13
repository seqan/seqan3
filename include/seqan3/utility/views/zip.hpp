// -----------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2022, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2022, Knut Reinert & MPI für molekulare Genetik
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at: https://github.com/seqan/seqan3/blob/master/LICENSE.md
// -----------------------------------------------------------------------------------------------------

/*!\file
 * \author Enrico Seiler <enrico.seiler AT fu-berlin.de>
 * \brief Provides seqan3::views::zip.
 */

#pragma once

#include <algorithm>
#include <functional>
#include <ranges>

#include <seqan3/core/detail/all_view.hpp>
#include <seqan3/core/range/detail/adaptor_from_functor.hpp>
#include <seqan3/utility/tuple/common_tuple.hpp>
#include <seqan3/utility/views/detail/helper.hpp>

//!\cond
// Contains helpers for seqan3::views::zip.
namespace seqan3::detail::zip
{

// https://eel.is/c++draft/range.zip#concept:zip-is-common
template <typename... range_ts>
concept zip_is_common = (sizeof...(range_ts) == 1 && (std::ranges::common_range<range_ts> && ...))
                     || (!(std::ranges::bidirectional_range<range_ts> && ...)
                         && (std::ranges::common_range<range_ts> && ...))
                     || ((std::ranges::random_access_range<range_ts> && ...)
                         && (std::ranges::sized_range<range_ts> && ...));

// Helper for tuple_or_pair
template <typename... ts>
struct tuple_or_pair_impl;

// Helper for tuple_or_pair
template <typename... ts>
    requires (sizeof...(ts) != 2)
struct tuple_or_pair_impl<ts...>
{
    using type = seqan3::common_tuple<ts...>;
};

// Helper for tuple_or_pair
template <typename... ts>
    requires (sizeof...(ts) == 2)
struct tuple_or_pair_impl<ts...>
{
    using type = seqan3::common_pair<ts...>;
};

// https://eel.is/c++draft/range.zip#view-1
template <typename... ts>
using tuple_or_pair = tuple_or_pair_impl<ts...>::type;

// std::abs has problems with ambiguous overloads
template <typename t>
constexpr t abs(t && v) noexcept
{
    if constexpr (std::is_signed_v<t>)
        return v < 0 ? -v : v;
    else
        return v;
}

// Returns a new tuple containing the result of applying a function to each tuple element.
// https://eel.is/c++draft/range.zip.view
template <typename fun_t, typename tuple_t>
constexpr auto tuple_transform(fun_t && f, tuple_t && tuple)
{
    return std::apply(
        [&]<typename... ts>(ts &&... elements)
        {
            return tuple_or_pair<std::invoke_result_t<fun_t &, ts>...>{std::invoke(f, std::forward<ts>(elements))...};
        },
        std::forward<tuple_t>(tuple));
}

// Applies a function to each tuple element.
// https://eel.is/c++draft/range.zip.view
template <typename fun_t, typename tuple_t>
constexpr void tuple_for_each(fun_t && f, tuple_t && tuple)
{
    std::apply(
        [&]<typename... ts>(ts &&... elements)
        {
            (std::invoke(f, std::forward<ts>(elements)), ...);
        },
        std::forward<tuple_t>(tuple));
}

template <bool is_const, typename... range_ts>
concept all_random_access = (std::ranges::random_access_range<view_helper::maybe_const<is_const, range_ts>> && ...);

template <bool is_const, typename... range_ts>
concept all_bidirectional = (std::ranges::bidirectional_range<view_helper::maybe_const<is_const, range_ts>> && ...);

template <bool is_const, typename... range_ts>
concept all_forward = (std::ranges::forward_range<view_helper::maybe_const<is_const, range_ts>> && ...);

// Pre cpp20-iterators: Infer iterator_category from modelled concepts.
#if defined(__GNUC__) && (__GNUC__ == 10)
template <bool is_const, typename... range_ts>
concept all_contiguous = (std::ranges::contiguous_range<view_helper::maybe_const<is_const, range_ts>> && ...);

template <bool Const, typename... Views>
struct iterator_category_t
{
    using iterator_category = std::conditional_t<
        all_contiguous<Const, Views...>,
        std::contiguous_iterator_tag,
        std::conditional_t<
            all_random_access<Const, Views...>,
            std::random_access_iterator_tag,
            std::conditional_t<
                all_bidirectional<Const, Views...>,
                std::bidirectional_iterator_tag,
                std::conditional_t<all_forward<Const, Views...>, std::forward_iterator_tag, std::input_iterator_tag>>>>;
};
#else // cpp20 iterators
template <bool>
struct iterator_category_t;
template <>
struct iterator_category_t<true>
{
    using iterator_category = std::input_iterator_tag;
};
template <>
struct iterator_category_t<false>
{};
#endif

} // namespace seqan3::detail::zip

namespace seqan3::detail
{

template <std::ranges::input_range... Views>
    requires (std::ranges::view<Views> && ...) && (sizeof...(Views) > 0)
class zip_view : public std::ranges::view_interface<zip_view<Views...>>
{
    seqan3::common_tuple<Views...> views_;

    template <bool>
    class iterator;

    template <bool>
    class sentinel;

public:
    zip_view()
        requires (std::is_default_constructible_v<Views> && ...)
    = default;
    constexpr explicit zip_view(Views... views) : views_(std::move(views)...)
    {}

    constexpr auto begin()
        requires (!(view_helper::simple_view<Views> && ...))
    {
        return iterator<false>(zip::tuple_transform(std::ranges::begin, views_));
    }

    constexpr auto begin() const
        requires (std::ranges::range<Views const> && ...)
    {
        return iterator<true>(zip::tuple_transform(std::ranges::begin, views_));
    }

    constexpr auto end()
        requires (!(view_helper::simple_view<Views> && ...))
    {
        if constexpr (!zip::zip_is_common<Views...>)
            return sentinel<false>(zip::tuple_transform(std::ranges::end, views_));
        else if constexpr ((std::ranges::random_access_range<Views> && ...))
            return begin() + std::iter_difference_t<iterator<false>>(size());
        else
            return iterator<false>(zip::tuple_transform(std::ranges::end, views_));
    }

    constexpr auto end() const
        requires (std::ranges::range<Views const> && ...)
    {
        if constexpr (!zip::zip_is_common<Views const...>)
            return sentinel<true>(zip::tuple_transform(std::ranges::end, views_));
        else if constexpr ((std::ranges::random_access_range<Views const> && ...))
            return begin() + std::iter_difference_t<iterator<true>>(size());
        else
            return iterator<true>(zip::tuple_transform(std::ranges::end, views_));
    }

    constexpr auto size()
        requires (std::ranges::sized_range<Views> && ...)
    {
        return std::apply(
            [](auto... sizes)
            {
                using common_size_t = std::make_unsigned_t<std::common_type_t<decltype(sizes)...>>;
                return std::ranges::min({static_cast<common_size_t>(sizes)...});
            },
            zip::tuple_transform(std::ranges::size, views_));
    }

    constexpr auto size() const
        requires (std::ranges::sized_range<Views const> && ...)
    {
        return std::apply(
            [](auto... sizes)
            {
                using common_size_t = std::make_unsigned_t<std::common_type_t<decltype(sizes)...>>;
                return std::ranges::min({static_cast<common_size_t>(sizes)...});
            },
            zip::tuple_transform(std::ranges::size, views_));
    }
};

template <typename... range_ts>
zip_view(range_ts &&...) -> zip_view<seqan3::detail::all_t<range_ts>...>;

template <std::ranges::input_range... Views>
    requires (std::ranges::view<Views> && ...) && (sizeof...(Views) > 0)
template <bool Const>
#if defined(__GNUC__) && (__GNUC__ == 10) // cpp17 iterators
class zip_view<Views...>::iterator : public zip::iterator_category_t<Const, Views...>
#else // cpp20 iterators
class zip_view<Views...>::iterator : public zip::iterator_category_t<zip::all_forward<Const, Views...>>
#endif
{
private:
    constexpr explicit iterator(
        zip::tuple_or_pair<std::ranges::iterator_t<view_helper::maybe_const<Const, Views>>...> current) :
        current_(std::move(current))
    {}

    friend class zip_view<Views...>;

public:
    // Exposition-only. Is public for access via friend operator== of the sentinel.
    zip::tuple_or_pair<std::ranges::iterator_t<view_helper::maybe_const<Const, Views>>...> current_;

    using iterator_concept = std::conditional_t<
        zip::all_random_access<Const, Views...>,
        std::random_access_iterator_tag,
        std::conditional_t<
            zip::all_bidirectional<Const, Views...>,
            std::bidirectional_iterator_tag,
            std::conditional_t<zip::all_forward<Const, Views...>, std::forward_iterator_tag, std::input_iterator_tag>>>;
    using value_type = zip::tuple_or_pair<std::ranges::range_value_t<view_helper::maybe_const<Const, Views>>...>;
    using difference_type =
        std::common_type_t<std::ranges::range_difference_t<view_helper::maybe_const<Const, Views>>...>;

    iterator() = default;
    constexpr iterator(iterator<!Const> i)
        requires Const
              && (std::convertible_to<std::ranges::iterator_t<Views>,
                                      std::ranges::iterator_t<view_helper::maybe_const<Const, Views>>>
                  && ...)
    : current_(std::move(i.current))
    {}

    constexpr auto operator*() const
    {
        return zip::tuple_transform(
            [](auto & i) -> decltype(auto)
            {
                return *i;
            },
            current_);
    }

    constexpr iterator & operator++()
    {
        zip::tuple_for_each(
            [](auto & i)
            {
                ++i;
            },
            current_);
        return *this;
    }

    constexpr void operator++(int)
    {
        ++*this;
    }

    constexpr iterator operator++(int)
        requires zip::all_forward<Const, Views...>
    {
        auto tmp = *this;
        ++*this;
        return tmp;
    }

    constexpr iterator & operator--()
        requires zip::all_bidirectional<Const, Views...>
    {
        zip::tuple_for_each(
            [](auto & i)
            {
                --i;
            },
            current_);
        return *this;
    }

    constexpr iterator operator--(int)
        requires zip::all_bidirectional<Const, Views...>
    {
        auto tmp = *this;
        --*this;
        return tmp;
    }

    constexpr iterator & operator+=(difference_type x)
        requires zip::all_random_access<Const, Views...>
    {
        zip::tuple_for_each(
            [&]<typename I>(I & i)
            {
                i += std::iter_difference_t<I>(x);
            },
            current_);
        return *this;
    }

    constexpr iterator & operator-=(difference_type x)
        requires zip::all_random_access<Const, Views...>
    {
        zip::tuple_for_each(
            [&]<typename I>(I & i)
            {
                i -= std::iter_difference_t<I>(x);
            },
            current_);
        return *this;
    }

    constexpr auto operator[](difference_type n) const
        requires zip::all_random_access<Const, Views...>
    {
        return zip::tuple_transform(
            [&]<typename I>(I & i) -> decltype(auto)
            {
                return i[std::iter_difference_t<I>(n)];
            },
            current_);
    }

    friend constexpr bool operator==(iterator const & x, iterator const & y)
        requires (std::equality_comparable<std::ranges::iterator_t<view_helper::maybe_const<Const, Views>>> && ...)
    {
        if constexpr (zip::all_bidirectional<Const, Views...>)
        {
            return x.current_ == y.current_;
        }
        else
        {
            return [&]<size_t... N>(std::integer_sequence<size_t, N...>)
            {
                return ((std::get<N>(x.current_) == std::get<N>(y.current_)) || ...);
            }
            (std::index_sequence_for<Views...>{});
        }
    }

    friend constexpr bool operator<(iterator const & x, iterator const & y)
        requires zip::all_random_access<Const, Views...>
    {
        return x.current_ < y.current_;
    }

    friend constexpr bool operator>(iterator const & x, iterator const & y)
        requires zip::all_random_access<Const, Views...>
    {
        return y < x;
    }

    friend constexpr bool operator<=(iterator const & x, iterator const & y)
        requires zip::all_random_access<Const, Views...>
    {
        return !(y < x);
    }

    friend constexpr bool operator>=(iterator const & x, iterator const & y)
        requires zip::all_random_access<Const, Views...>
    {
        return !(x < y);
    }

#ifdef __cpp_lib_three_way_comparison
    friend constexpr auto operator<=>(iterator const & x, iterator const & y)
        requires zip::all_random_access<Const, Views...>
              && (std::three_way_comparable<std::ranges::iterator_t<view_helper::maybe_const<Const, Views>>> && ...)
    {
        return x.current_ <=> y.current_;
    }
#endif

    friend constexpr iterator operator+(iterator const & i, difference_type n)
        requires zip::all_random_access<Const, Views...>
    {
        auto r = i;
        r += n;
        return r;
    }

    friend constexpr iterator operator+(difference_type n, iterator const & i)
        requires zip::all_random_access<Const, Views...>
    {
        return i + n;
    }

    friend constexpr iterator operator-(iterator const & i, difference_type n)
        requires zip::all_random_access<Const, Views...>
    {
        auto r = i;
        r -= n;
        return r;
    }

    friend constexpr difference_type operator-(iterator const & x, iterator const & y)
        requires (std::sized_sentinel_for<std::ranges::iterator_t<view_helper::maybe_const<Const, Views>>,
                                          std::ranges::iterator_t<view_helper::maybe_const<Const, Views>>>
                  && ...)
    {
        return [&]<size_t... N>(std::integer_sequence<size_t, N...>)
        {
            return std::ranges::min(
                {static_cast<difference_type>(std::get<N>(x.current_) - std::get<N>(y.current_))...},
                [](difference_type a, difference_type b)
                {
                    return zip::abs(b) < zip::abs(a);
                });
        }
        (std::index_sequence_for<Views...>{});
    }

    friend constexpr auto iter_move(iterator const & i) noexcept(
        (noexcept(std::ranges::iter_move(
             std::declval<std::ranges::iterator_t<view_helper::maybe_const<Const, Views>> const &>()))
         && ...)
        && (std::is_nothrow_move_constructible_v<
                std::ranges::range_rvalue_reference_t<view_helper::maybe_const<Const, Views>>>
            && ...))
    {
        return zip::tuple_transform(std::ranges::iter_move, i.current_);
    }

    friend constexpr void iter_swap(iterator const & l, iterator const & r) noexcept(
        (noexcept(std::ranges::iter_swap(
             std::declval<std::ranges::iterator_t<view_helper::maybe_const<Const, Views>> const &>(),
             std::declval<std::ranges::iterator_t<view_helper::maybe_const<Const, Views>> const &>()))
         && ...))
        requires (std::indirectly_swappable<std::ranges::iterator_t<view_helper::maybe_const<Const, Views>>> && ...)
    {
        [&]<size_t... N>(std::integer_sequence<size_t, N...>)
        {
            (std::ranges::iter_swap(std::get<N>(l.current_), std::get<N>(r.current)), ...);
        }
        (std::index_sequence_for<Views...>{});
    }
};

template <std::ranges::input_range... Views>
    requires (std::ranges::view<Views> && ...) && (sizeof...(Views) > 0)
template <bool Const>
class zip_view<Views...>::sentinel
{
private:
    constexpr explicit sentinel(
        zip::tuple_or_pair<std::ranges::sentinel_t<view_helper::maybe_const<Const, Views>>...> end) :
        end_(std::move(end))
    {}

    friend class zip_view<Views...>;

public:
    // Exposition-only. Is public such that it can be accessed by friends.
    zip::tuple_or_pair<std::ranges::sentinel_t<view_helper::maybe_const<Const, Views>>...> end_;

    sentinel() = default;
    constexpr sentinel(sentinel<!Const> i)
        requires Const
              && (std::convertible_to<std::ranges::sentinel_t<Views>,
                                      std::ranges::sentinel_t<view_helper::maybe_const<Const, Views>>>
                  && ...)
    : end_(std::move(i.end_))
    {}

    template <bool OtherConst>
        requires (std::sentinel_for<std::ranges::sentinel_t<view_helper::maybe_const<Const, Views>>,
                                    std::ranges::iterator_t<view_helper::maybe_const<OtherConst, Views>>>
                  && ...)
    friend constexpr bool operator==(iterator<OtherConst> const & x, sentinel const & y)
    {
        return [&]<size_t... N>(std::integer_sequence<size_t, N...>)
        {
            return ((std::get<N>(x.current_) == std::get<N>(y.end_)) || ...);
        }
        (std::index_sequence_for<Views...>{});
    }

    template <bool OtherConst>
        requires (std::sized_sentinel_for<std::ranges::sentinel_t<view_helper::maybe_const<Const, Views>>,
                                          std::ranges::iterator_t<view_helper::maybe_const<OtherConst, Views>>>
                  && ...)
    friend constexpr std::common_type_t<std::ranges::range_difference_t<view_helper::maybe_const<OtherConst, Views>>...>
    operator-(iterator<OtherConst> const & x, sentinel const & y)
    {
        using return_t =
            std::common_type_t<std::ranges::range_difference_t<view_helper::maybe_const<OtherConst, Views>>...>;
        return [&]<size_t... N>(std::integer_sequence<size_t, N...>)
        {
            return std::ranges::min({static_cast<return_t>(std::get<N>(x.current_) - std::get<N>(y.end_))...},
                                    [](return_t a, return_t b)
                                    {
                                        return zip::abs(b) < zip::abs(a);
                                    });
        }
        (std::index_sequence_for<Views...>{});
    }

    template <bool OtherConst>
        requires (std::sized_sentinel_for<std::ranges::sentinel_t<view_helper::maybe_const<Const, Views>>,
                                          std::ranges::iterator_t<view_helper::maybe_const<OtherConst, Views>>>
                  && ...)
    friend constexpr std::common_type_t<std::ranges::range_difference_t<view_helper::maybe_const<OtherConst, Views>>...>
    operator-(sentinel const & y, iterator<OtherConst> const & x)
    {
        return -(x - y);
    }
};

struct zip_fn
{
    template <typename urng_t>
    constexpr auto operator()(urng_t && rng) const
    {
        return adaptor_from_functor{*this, std::forward<urng_t>(rng)};
    }

    template <typename... urng_ts>
        requires (sizeof...(urng_ts) == 0)
    constexpr auto operator()(urng_ts &&... ranges) const
    {
        return std::views::empty<seqan3::common_tuple<>>;
    }

    template <typename... urng_ts>
        requires (sizeof...(urng_ts) > 1)
    constexpr auto operator()(urng_ts &&... ranges) const
    {
        return zip_view{std::forward<urng_ts>(ranges)...};
    }
};

} // namespace seqan3::detail
//!\endcond

namespace seqan3::views
{

/*!\brief A view adaptor that takes several views and returns tuple-like values from every i-th element of each view.
 * \ingroup utility_views
 * \noapi{This is a implementation of the C++23 zip_view. It will be replaced with std::views::zip.}
 * \sa https://en.cppreference.com/w/cpp/ranges/zip_view
 */
inline constexpr auto zip = seqan3::detail::zip_fn{};

} // namespace seqan3::views
