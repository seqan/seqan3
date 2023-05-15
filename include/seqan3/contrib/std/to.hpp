// -----------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2023, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2023, Knut Reinert & MPI für molekulare Genetik
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at: https://github.com/seqan/seqan-std/blob/main/LICENSE
// -----------------------------------------------------------------------------------------------------

/*!\file
 * \author Enrico Seiler <enrico.seiler AT fu-berlin.de>
 * \brief Provides seqan::std::ranges::to.
 */

// File might be included from multiple libraries.
#ifndef SEQAN_STD_TO
#define SEQAN_STD_TO

#include <ranges>

#ifdef __cpp_lib_ranges_to_container

namespace seqan::std
{

using ::std::from_range_t;

inline constexpr from_range_t from_range{};

} // namespace seqan::std

namespace seqan::std::ranges
{

using ::std::ranges::to;

} // namespace seqan::std::ranges

#else

#    include <algorithm>

#    include "detail/adaptor_from_functor.hpp"

namespace seqan::std
{

struct from_range_t
{
    explicit from_range_t() = default;
};

inline constexpr from_range_t from_range{};

} // namespace seqan::std

namespace seqan::std::detail::to
{

// clang-format off
template <class Container>
constexpr bool reservable_container =
    ::std::ranges::sized_range<Container> &&
    requires (Container & c, ::std::ranges::range_size_t<Container> n)
    {
        c.reserve(n);
        { c.capacity() } -> ::std::same_as<decltype(n)>;
        { c.max_size() } -> ::std::same_as<decltype(n)>;
    };

template <class Container, class Ref>
constexpr bool container_insertable =
    requires (Container & c, Ref && ref)
    {
        requires
        (
            requires { c.push_back(::std::forward<Ref>(ref)); } ||
            requires { c.insert(c.end(), ::std::forward<Ref>(ref)); }
        );
    };
// clang-format on

template <class Ref, class Container>
constexpr auto container_inserter(Container & c)
{
    if constexpr (requires { c.push_back(::std::declval<Ref>()); })
        return ::std::back_inserter(c);
    else
        return ::std::inserter(c, c.end());
}

template <::std::ranges::input_range R>
struct input_iterator
{
    using iterator_category = ::std::input_iterator_tag;
    using value_type = ::std::ranges::range_value_t<R>;
    using difference_type = ::std::ptrdiff_t;
    using pointer = ::std::add_pointer_t<::std::ranges::range_reference_t<R>>;
    using reference = ::std::ranges::range_reference_t<R>;

    reference operator*() const;
    pointer operator->() const;
    input_iterator & operator++();
    input_iterator operator++(int);
    bool operator==(input_iterator const &) const;
};

template <typename It>
concept has_input_iterator_category =
    ::std::derived_from<typename ::std::iterator_traits<It>::iterator_category, ::std::input_iterator_tag>;

} // namespace seqan::std::detail::to

namespace seqan::std::ranges
{

template <class C, ::std::ranges::input_range R, class... Args>
    requires (!::std::ranges::view<C>)
constexpr C to(R && r, Args &&... args)
{
    if constexpr (!::std::ranges::input_range<C>
                  || ::std::convertible_to<::std::ranges::range_reference_t<R>, ::std::ranges::range_value_t<C>>)
    {
        if constexpr (::std::constructible_from<C, R, Args...>)
            return C(::std::forward<R>(r), ::std::forward<Args>(args)...);
        else if constexpr (::std::constructible_from<C, seqan::std::from_range_t, R, Args...>)
            return C(seqan::std::from_range, ::std::forward<R>(r), ::std::forward<Args>(args)...);
        else if constexpr (::std::ranges::common_range<R>
                           && seqan::std::detail::to::has_input_iterator_category<::std::ranges::iterator_t<R>>
                           && ::std::constructible_from<C,
                                                        ::std::ranges::iterator_t<R>,
                                                        ::std::ranges::sentinel_t<R>,
                                                        Args...>)
            return C(::std::ranges::begin(r), ::std::ranges::end(r), ::std::forward<Args>(args)...);
        else if constexpr (::std::constructible_from<C, Args...>
                           && seqan::std::detail::to::container_insertable<C, ::std::ranges::range_reference_t<R>>)
        {
            C c(::std::forward<Args>(args)...);
            if constexpr (::std::ranges::sized_range<R> && seqan::std::detail::to::reservable_container<C>)
                c.reserve(static_cast<::std::ranges::range_size_t<C>>(::std::ranges::size(r)));
            ::std::ranges::copy(r, seqan::std::detail::to::container_inserter<::std::ranges::range_reference_t<R>>(c));
            return c;
        }
    }
    else if constexpr (::std::ranges::input_range<::std::ranges::range_reference_t<R>>)
        return to<C>(r
                         | ::std::views::transform(
                             [](auto && elem)
                             {
                                 return to<::std::ranges::range_value_t<C>>(::std::forward<decltype(elem)>(elem));
                             }),
                     ::std::forward<Args>(args)...);
    else
        __builtin_unreachable();
}

template <template <class...> class C, ::std::ranges::input_range R, class... Args>
constexpr auto to(R && r, Args &&... args)
{
    if constexpr (requires { C(::std::declval<R>(), ::std::declval<Args>()...); })
        return to<decltype(C(::std::declval<R>(), ::std::declval<Args>()...))>(::std::forward<R>(r),
                                                                               ::std::forward<Args>(args)...);
    else if constexpr (requires { C(seqan::std::from_range, ::std::declval<R>(), ::std::declval<Args>()...); })
        return to<decltype(C(seqan::std::from_range, ::std::declval<R>(), ::std::declval<Args>()...))>(
            ::std::forward<R>(r),
            ::std::forward<Args>(args)...);
    else if constexpr (requires {
                           C(::std::declval<seqan::std::detail::to::input_iterator<R>>(),
                             ::std::declval<seqan::std::detail::to::input_iterator<R>>(),
                             ::std::declval<Args>()...);
                       })
        return to<decltype(C(::std::declval<seqan::std::detail::to::input_iterator<R>>(),
                             ::std::declval<seqan::std::detail::to::input_iterator<R>>(),
                             ::std::declval<Args>()...))>(::std::forward<R>(r), ::std::forward<Args>(args)...);
    else
        __builtin_unreachable();
}

} // namespace seqan::std::ranges

namespace seqan::std::detail::to
{

template <class C>
    requires (!::std::ranges::view<C>)
struct to_fn1
{
    template <class... Args>
    constexpr auto operator()(Args &&... args) const
    {
        return seqan::std::detail::adaptor_from_functor{*this, ::std::forward<Args>(args)...};
    }

    template <::std::ranges::input_range R, class... Args>
    constexpr auto operator()(R && r, Args &&... args) const
    {
        return seqan::std::ranges::to<C>(::std::forward<R>(r), ::std::forward<Args>(args)...);
    }
};

template <template <class...> class C>
struct to_fn2
{
    template <class... Args>
    constexpr auto operator()(Args &&... args) const
    {
        return seqan::std::detail::adaptor_from_functor{*this, ::std::forward<Args>(args)...};
    }

    template <::std::ranges::input_range R, class... Args>
    constexpr auto operator()(R && r, Args &&... args) const
    {
        return seqan::std::ranges::to<C>(::std::forward<R>(r), ::std::forward<Args>(args)...);
    }
};

} // namespace seqan::std::detail::to

namespace seqan::std::ranges
{

template <class C, class... Args>
    requires (!::std::ranges::view<C>)
constexpr auto to(Args &&... args)
{
    return seqan::std::detail::to::to_fn1<C>{}(::std::forward<Args>(args)...);
}

template <template <class...> class C, class... Args>
constexpr auto to(Args &&... args)
{
    return seqan::std::detail::to::to_fn2<C>{}(::std::forward<Args>(args)...);
}

} // namespace seqan::std::ranges

#endif // ifdef __cpp_lib_ranges_to_container

#endif // SEQAN_STD_TO
