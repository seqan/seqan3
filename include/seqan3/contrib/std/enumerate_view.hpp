// SPDX-FileCopyrightText: 2006-2025 Knut Reinert & Freie Universität Berlin
// SPDX-FileCopyrightText: 2016-2025 Knut Reinert & MPI für molekulare Genetik
// SPDX-License-Identifier: BSD-3-Clause

/*!\file
 * \author Enrico Seiler <enrico.seiler AT fu-berlin.de>
 * \brief Provides seqan::stl::views::enumerate.
 */

// File might be included from multiple libraries.
#ifndef SEQAN_STD_ENUMERATE_VIEW
#define SEQAN_STD_ENUMERATE_VIEW

#include <cassert>
#include <ranges>

#ifdef __cpp_lib_ranges_enumerate

namespace seqan::stl::views
{

using std::ranges::views::enumerate;

} // namespace seqan::stl::views

#else

#    include "all_view.hpp"
#    include "concepts.hpp"
#    include "detail/adaptor_from_functor.hpp"
#    include "detail/compiler_definitions.hpp"
#    include "detail/exposition_only.hpp"

namespace seqan::stl::ranges
{

template <std::ranges::view V>
    requires seqan::stl::detail::range_with_movable_references<V>
class enumerate_view : public std::ranges::view_interface<enumerate_view<V>>
{
    // clang-format off
SEQAN_STD_NESTED_VISIBILITY
    // clang-format on
    V base_ = V();

    template <bool Const>
    class iterator;

    template <bool Const>
    class sentinel;

public:
    constexpr enumerate_view()
        requires std::default_initializable<V>
    = default;

    constexpr explicit enumerate_view(V base) : base_{std::move(base)}
    {}

    constexpr auto begin()
        requires (!seqan::stl::detail::simple_view<V>)
    {
        return iterator<false>{std::ranges::begin(base_), 0};
    }
    constexpr auto begin() const
        requires seqan::stl::detail::range_with_movable_references<V const>
    {
        return iterator<true>{std::ranges::begin(base_), 0};
    }

    constexpr auto end()
        requires (!seqan::stl::detail::simple_view<V>)
    {
        if constexpr (std::ranges::forward_range<V> && std::ranges::common_range<V> && std::ranges::sized_range<V>)
            return iterator<false>{std::ranges::end(base_), std::ranges::distance(base_)};
        else
            return sentinel<false>{std::ranges::end(base_)};
    }
    constexpr auto end() const
        requires seqan::stl::detail::range_with_movable_references<V const>
    {
        if constexpr (std::ranges::forward_range<V const> && std::ranges::common_range<V const>
                      && std::ranges::sized_range<V const>)
            return iterator<true>{std::ranges::end(base_), std::ranges::distance(base_)};
        else
            return sentinel<true>{std::ranges::end(base_)};
    }

    constexpr auto size()
        requires std::ranges::sized_range<V>
    {
        return std::ranges::size(base_);
    }
    constexpr auto size() const
        requires std::ranges::sized_range<V const>
    {
        return std::ranges::size(base_);
    }

    constexpr V base() const &
        requires std::copy_constructible<V>
    {
        return base_;
    }
    constexpr V base() &&
    {
        return std::move(base_);
    }
};

template <class R>
enumerate_view(R &&) -> enumerate_view<seqan::stl::views::all_t<R>>;

template <std::ranges::view V>
    requires seqan::stl::detail::range_with_movable_references<V>
template <bool Const>
class enumerate_view<V>::iterator
{
    friend enumerate_view;
    using Base = seqan::stl::detail::maybe_const<Const, V>;

public:
    using iterator_category = std::input_iterator_tag;
    using iterator_concept = std::conditional_t<
        std::ranges::random_access_range<Base>,
        std::random_access_iterator_tag,
        std::conditional_t<
            std::ranges::bidirectional_range<Base>,
            std::bidirectional_iterator_tag,
            std::conditional_t<std::ranges::forward_range<Base>, std::forward_iterator_tag, std::input_iterator_tag>>>;
    using difference_type = std::ranges::range_difference_t<Base>;
    using value_type = std::tuple<std::ranges::range_difference_t<Base>, std::ranges::range_value_t<Base>>;

private:
    using reference_type = std::tuple<std::ranges::range_difference_t<Base>, std::ranges::range_reference_t<Base>>;
    std::ranges::iterator_t<Base> current_ = std::ranges::iterator_t<Base>{};
    difference_type pos_ = 0;

    constexpr explicit iterator(std::ranges::iterator_t<Base> current, difference_type pos) :
        current_{std::move(current)},
        pos_{pos}
    {}

public:
    iterator()
        requires std::default_initializable<std::ranges::iterator_t<Base>>
    = default;
    constexpr iterator(iterator<!Const> i)
        requires Const && std::convertible_to<std::ranges::iterator_t<V>, std::ranges::iterator_t<Base>>
        : current_{std::move(i.current_)}, pos_{i.pos_}
    {}

    constexpr std::ranges::iterator_t<Base> const & base() const & noexcept
    {
        return current_;
    }
    constexpr std::ranges::iterator_t<Base> base() &&
    {
        return std::move(current_);
    }

    constexpr difference_type index() const noexcept
    {
        return pos_;
    }

    constexpr auto operator*() const
    {
        return reference_type{pos_, *current_};
    }

    constexpr iterator & operator++()
    {
        ++current_;
        ++pos_;
        return *this;
    }
    constexpr void operator++(int)
    {
        ++*this;
    }
    constexpr iterator operator++(int)
        requires std::ranges::forward_range<Base>
    {
        auto temp = *this;
        ++*this;
        return temp;
    }

    constexpr iterator & operator--()
        requires std::ranges::bidirectional_range<Base>
    {
        --current_;
        --pos_;
        return *this;
    }
    constexpr iterator operator--(int)
        requires std::ranges::bidirectional_range<Base>
    {
        auto temp = *this;
        --*this;
        return temp;
    }

    constexpr iterator & operator+=(difference_type n)
        requires std::ranges::random_access_range<Base>
    {
        current_ += n;
        pos_ += n;
        return *this;
    }
    constexpr iterator & operator-=(difference_type n)
        requires std::ranges::random_access_range<Base>
    {
        current_ -= n;
        pos_ -= n;
        return *this;
    }

    constexpr auto operator[](difference_type x) const
        requires std::ranges::random_access_range<Base>
    {
        return reference_type{pos_ + x, current_[x]};
    }

    friend constexpr bool operator==(iterator const & x, iterator const & y) noexcept
    {
        return x.pos_ == y.pos_;
    }
#    ifdef __cpp_lib_three_way_comparison
    friend constexpr std::strong_ordering operator<=>(iterator const & x, iterator const & y) noexcept
    {
        return x.pos_ <=> y.pos_;
    }
#    endif

    friend constexpr iterator operator+(iterator const & x, difference_type y)
        requires std::ranges::random_access_range<Base>
    {
        auto temp = x;
        temp += y;
        return temp;
    }
    friend constexpr iterator operator+(difference_type x, iterator const & y)
        requires std::ranges::random_access_range<Base>
    {
        return y + x;
    }
    friend constexpr iterator operator-(iterator const & x, difference_type y)
        requires std::ranges::random_access_range<Base>
    {
        auto temp = x;
        temp -= y;
        return temp;
    }
    friend constexpr difference_type operator-(iterator const & x, iterator const & y) noexcept
    {
        return x.pos_ - y.pos_;
    }

    friend constexpr auto
    iter_move(iterator const & i) noexcept(noexcept(std::ranges::iter_move(i.current_))
                                           && std::is_nothrow_move_constructible_v<std::ranges::range_value_t<Base>>)
    {
        return std::tuple<difference_type, std::ranges::range_rvalue_reference_t<Base>>{
            i.pos_,
            std::ranges::iter_move(i.current_)};
    }
};

template <std::ranges::view V>
    requires seqan::stl::detail::range_with_movable_references<V>
template <bool Const>
class enumerate_view<V>::sentinel
{
    friend enumerate_view;
    using Base = seqan::stl::detail::maybe_const<Const, V>;
    std::ranges::sentinel_t<Base> end_ = std::ranges::sentinel_t<Base>{};
    constexpr explicit sentinel(std::ranges::sentinel_t<Base> end) : end_{std::move(end)}
    {}

public:
    sentinel() = default;
    constexpr sentinel(sentinel<!Const> other)
        requires Const && std::convertible_to<std::ranges::sentinel_t<V>, std::ranges::sentinel_t<Base>>
        : end_{std::move(other.end_)}
    {}

    constexpr std::ranges::sentinel_t<Base> base() const
    {
        return end_;
    }

    template <bool OtherConst>
        requires std::sentinel_for<std::ranges::sentinel_t<Base>,
                                   std::ranges::iterator_t<seqan::stl::detail::maybe_const<OtherConst, V>>>
    friend constexpr bool operator==(iterator<OtherConst> const & x, sentinel const & y) noexcept
    {
        return x.current_ == y.end_;
    }

    template <bool OtherConst>
        requires std::sentinel_for<std::ranges::sentinel_t<Base>,
                                   std::ranges::iterator_t<seqan::stl::detail::maybe_const<OtherConst, V>>>
    friend constexpr std::ranges::range_difference_t<seqan::stl::detail::maybe_const<OtherConst, V>>
    operator-(iterator<OtherConst> const & x, sentinel const & y)
    {
        return x.current_ - y.end_;
    }

    template <bool OtherConst>
        requires std::sized_sentinel_for<std::ranges::sentinel_t<Base>,
                                         std::ranges::iterator_t<seqan::stl::detail::maybe_const<OtherConst, V>>>
    friend constexpr std::ranges::range_difference_t<seqan::stl::detail::maybe_const<OtherConst, V>>
    operator-(sentinel const & x, iterator<OtherConst> const & y)
    {
        return x.end_ - y.current_;
    }
};

struct enumerate_fn
{
    template <seqan::stl::ranges::viewable_range Range>
    constexpr auto operator()(Range && range) const
    {
        return enumerate_view{std::forward<Range>(range)};
    }
};

} // namespace seqan::stl::ranges

namespace seqan::stl::views
{

inline constexpr auto enumerate = seqan::stl::ranges::enumerate_fn{};

} // namespace seqan::stl::views

#endif // ifdef __cpp_lib_ranges_enumerate

#endif // SEQAN_STD_ENUMERATE_VIEW
