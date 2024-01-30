// SPDX-FileCopyrightText: 2006-2024 Knut Reinert & Freie Universität Berlin
// SPDX-FileCopyrightText: 2016-2024 Knut Reinert & MPI für molekulare Genetik
// SPDX-License-Identifier: BSD-3-Clause

/*!\file
 * \author Enrico Seiler <enrico.seiler AT fu-berlin.de>
 * \brief Provides seqan::stl::views::chunk.
 */

// File might be included from multiple libraries.
#ifndef SEQAN_STD_CHUNK_VIEW
#define SEQAN_STD_CHUNK_VIEW

#include <cassert>
#include <ranges>

#ifdef __cpp_lib_ranges_chunk

namespace seqan::stl::views
{

using std::ranges::views::chunk;

} // namespace seqan::stl::views

#else

#    include <algorithm>

#    include "all_view.hpp"
#    include "concepts.hpp"
#    include "detail/adaptor_from_functor.hpp"
#    include "detail/compiler_definitions.hpp"
#    include "detail/exposition_only.hpp"
#    include "detail/non_propagating_cache.hpp"

namespace seqan::stl::detail::chunk
{

template <class I>
constexpr I div_ceil(I num, I denom)
{
    I r = num / denom;
    if (num % denom)
        ++r;
    return r;
}

template <std::integral T>
constexpr auto to_unsigned_like(T v) noexcept
{
    return static_cast<std::make_unsigned_t<T>>(v);
}

} // namespace seqan::stl::detail::chunk

namespace seqan::stl::ranges
{

template <std::ranges::view V>
class chunk_view;

template <std::ranges::view V>
    requires std::ranges::input_range<V>
class chunk_view<V> : public std::ranges::view_interface<chunk_view<V>>
{
    // clang-format off
SEQAN_STD_NESTED_VISIBILITY
    // clang-format on
    V base_;
    std::ranges::range_difference_t<V> n_;
    std::ranges::range_difference_t<V> remainder_ = 0;

    seqan::stl::detail::non_propagating_cache<std::ranges::iterator_t<V>> current_;

private:
    class outer_iterator;

    class inner_iterator;

public:
    chunk_view()
        requires std::default_initializable<V>
    = default;

    constexpr explicit chunk_view(V base, std::ranges::range_difference_t<V> n) : base_{std::move(base)}, n_{n}
    {
        assert(n > 0);
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

    constexpr outer_iterator begin()
    {
        current_ = std::ranges::begin(base_);
        remainder_ = n_;
        return outer_iterator{*this};
    }
    constexpr std::default_sentinel_t end() const noexcept
    {
        return std::default_sentinel;
    }

    constexpr auto size()
        requires std::ranges::sized_range<V>
    {
        return seqan::stl::detail::chunk::to_unsigned_like(
            seqan::stl::detail::chunk::div_ceil(std::ranges::distance(base_), n_));
    }
    constexpr auto size() const
        requires std::ranges::sized_range<V const>
    {
        return seqan::stl::detail::chunk::to_unsigned_like(
            seqan::stl::detail::chunk::div_ceil(std::ranges::distance(base_), n_));
    }
};

template <class R>
chunk_view(R &&, std::ranges::range_difference_t<R>) -> chunk_view<seqan::stl::views::all_t<R>>;

template <std::ranges::view V>
    requires std::ranges::input_range<V>
class chunk_view<V>::outer_iterator
{
private:
    chunk_view * parent_;

    constexpr explicit outer_iterator(chunk_view & parent) : parent_{std::addressof(parent)}
    {}

    friend chunk_view;

public:
    using iterator_concept = std::input_iterator_tag;
    using difference_type = std::ranges::range_difference_t<V>;

    struct value_type;

    outer_iterator(outer_iterator &&) = default;
    outer_iterator & operator=(outer_iterator &&) = default;

    constexpr value_type operator*() const
    {
        assert(*this != std::default_sentinel);
        return value_type(*parent_);
    }
    constexpr outer_iterator & operator++()
    {
        assert(*this != std::default_sentinel);
        std::ranges::advance(*parent_->current_, parent_->remainder_, std::ranges::end(parent_->base_));
        parent_->remainder_ = parent_->n_;
        return *this;
    }
    constexpr void operator++(int)
    {
        ++*this;
    }

    friend constexpr bool operator==(outer_iterator const & x, std::default_sentinel_t)
    {
        return *x.parent_->current_ == std::ranges::end(x.parent_->base_) && x.parent_->remainder_ != 0;
    }

    friend constexpr difference_type operator-(std::default_sentinel_t, outer_iterator const & x)
        requires std::sized_sentinel_for<std::ranges::sentinel_t<V>, std::ranges::iterator_t<V>>
    {
        auto const dist = std::ranges::end(x.parent_->base_) - *x.parent_->current_;
        if (dist < x.parent_->remainder_)
        {
            return dist == 0 ? 0 : 1;
        }
        return seqan::stl::detail::chunk::div_ceil(dist - x.parent_->remainder_, x.parent_->n_) + 1;
    }
    friend constexpr difference_type operator-(outer_iterator const & x, std::default_sentinel_t y)
        requires std::sized_sentinel_for<std::ranges::sentinel_t<V>, std::ranges::iterator_t<V>>
    {
        return -(y - x);
    }
};

template <std::ranges::view V>
    requires std::ranges::input_range<V>
struct chunk_view<V>::outer_iterator::value_type : std::ranges::view_interface<value_type>
{
private:
    chunk_view * parent_;

    constexpr explicit value_type(chunk_view & parent) noexcept : parent_{std::addressof(parent)}
    {}

    friend chunk_view;

public:
    constexpr inner_iterator begin() const noexcept
    {
        return inner_iterator{*parent_};
    }
    constexpr std::default_sentinel_t end() const noexcept
    {
        return std::default_sentinel;
    }

    constexpr auto size() const
        requires std::sized_sentinel_for<std::ranges::sentinel_t<V>, std::ranges::iterator_t<V>>
    {
        return seqan::stl::detail::chunk::to_unsigned_like(
            std::ranges::min(parent_->remainder_, std::ranges::end(parent_->base_) - *parent_->current_));
    }
};

template <std::ranges::view V>
    requires std::ranges::input_range<V>
class chunk_view<V>::inner_iterator
{
private:
    chunk_view * parent_;

    constexpr explicit inner_iterator(chunk_view & parent) noexcept : parent_{std::addressof(parent)}
    {}

    friend chunk_view;

public:
    using iterator_concept = std::input_iterator_tag;
    using difference_type = std::ranges::range_difference_t<V>;
    using value_type = std::ranges::range_value_t<V>;

    inner_iterator(inner_iterator &&) = default;
    inner_iterator & operator=(inner_iterator &&) = default;

    constexpr std::ranges::iterator_t<V> const & base() const &
    {
        return *parent_->current_;
    }

    constexpr std::ranges::range_reference_t<V> operator*() const
    {
        assert(*this != std::default_sentinel);
        return **parent_->current_;
    }
    constexpr inner_iterator & operator++()
    {
        assert(*this != std::default_sentinel);
        ++*parent_->current_;
        if (*parent_->current_ == std::ranges::end(parent_->base_))
            parent_->remainder_ = 0;
        else
            --parent_->remainder_;
        return *this;
    }
    constexpr void operator++(int)
    {
        return ++*this;
    }

    friend constexpr bool operator==(inner_iterator const & x, std::default_sentinel_t)
    {
        return x.parent_->remainder_ == 0;
    }

    friend constexpr difference_type operator-(std::default_sentinel_t, inner_iterator const & x)
        requires std::sized_sentinel_for<std::ranges::sentinel_t<V>, std::ranges::iterator_t<V>>
    {
        std::ranges::min(x.parent_->remainder_, std::ranges::end(x.parent_->base_) - *x.parent_->current_);
    }
    friend constexpr difference_type operator-(inner_iterator const & x, std::default_sentinel_t y)
        requires std::sized_sentinel_for<std::ranges::sentinel_t<V>, std::ranges::iterator_t<V>>
    {
        return -(y - x);
    }

    friend constexpr std::ranges::range_rvalue_reference_t<V>
    iter_move(inner_iterator const & i) noexcept(noexcept(std::ranges::iter_move(*i.parent_->current_)))
    {
        return std::ranges::iter_move(*i.parent_->current_);
    }

    friend constexpr void iter_swap(inner_iterator const & x, inner_iterator const & y) noexcept(
        noexcept(std::ranges::iter_swap(*x.parent_->current_, *y.parent_->current_)))
        requires std::indirectly_swappable<std::ranges::iterator_t<V>>
    {
        std::ranges::iter_swap(*x.parent_->current_, *y.parent_->current_);
    }
};

template <std::ranges::view V>
    requires std::ranges::forward_range<V>
class chunk_view<V> : public std::ranges::view_interface<chunk_view<V>>
{
private:
    V base_;
    std::ranges::range_difference_t<V> n_;

    template <bool>
    class iterator;

public:
    chunk_view()
        requires std::default_initializable<V>
    = default;

    constexpr explicit chunk_view(V base, std::ranges::range_difference_t<V> n) : base_{std::move(base)}, n_{n}
    {
        assert(n > 0);
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

    constexpr auto begin()
        requires (!seqan::stl::detail::simple_view<V>)
    {
        return iterator<false>{this, std::ranges::begin(base_)};
    }

    constexpr auto begin() const
        requires std::ranges::forward_range<V const>
    {
        return iterator<true>{this, std::ranges::begin(base_)};
    }

    constexpr auto end()
        requires (!seqan::stl::detail::simple_view<V>)
    {
        if constexpr (std::ranges::common_range<V> && std::ranges::sized_range<V>)
        {
            auto missing = (n_ - std::ranges::distance(base_) % n_) % n_;
            return iterator<false>{this, std::ranges::end(base_), missing};
        }
        else if constexpr (std::ranges::common_range<V> && !std::ranges::bidirectional_range<V>)
        {
            return iterator<false>{this, std::ranges::end(base_)};
        }
        else
        {
            return std::default_sentinel;
        }
    }

    constexpr auto end() const
        requires std::ranges::forward_range<V const>
    {
        if constexpr (std::ranges::common_range<V const> && std::ranges::sized_range<V const>)
        {
            auto missing = (n_ - std::ranges::distance(base_) % n_) % n_;
            return iterator<true>{this, std::ranges::end(base_), missing};
        }
        else if constexpr (std::ranges::common_range<V const> && !std::ranges::bidirectional_range<V const>)
        {
            return iterator<true>{this, std::ranges::end(base_)};
        }
        else
        {
            return std::default_sentinel;
        }
    }

    constexpr auto size()
        requires std::ranges::sized_range<V>
    {
        return seqan::stl::detail::chunk::to_unsigned_like(
            seqan::stl::detail::chunk::div_ceil(std::ranges::distance(base_), n_));
    }
    constexpr auto size() const
        requires std::ranges::sized_range<V const>
    {
        return seqan::stl::detail::chunk::to_unsigned_like(
            seqan::stl::detail::chunk::div_ceil(std::ranges::distance(base_), n_));
    }
};

template <std::ranges::view V>
    requires std::ranges::forward_range<V>
template <bool Const>
class chunk_view<V>::iterator
{
private:
    using Parent = seqan::stl::detail::maybe_const<Const, chunk_view>;
    using Base = seqan::stl::detail::maybe_const<Const, V>;

    std::ranges::iterator_t<Base> current_ = std::ranges::iterator_t<Base>{};
    std::ranges::sentinel_t<Base> end_ = std::ranges::sentinel_t<Base>{};
    std::ranges::range_difference_t<Base> n_ = 0;
    std::ranges::range_difference_t<Base> missing_ = 0;

    constexpr iterator(Parent * parent,
                       std::ranges::iterator_t<Base> current,
                       std::ranges::range_difference_t<Base> missing = 0) :
        current_{current},
        end_{std::ranges::end(parent->base_)},
        n_{parent->n_},
        missing_{missing}
    {}

    friend chunk_view;

public:
    using iterator_category = std::input_iterator_tag;
    using iterator_concept = std::conditional_t<std::ranges::random_access_range<Base>,
                                                std::random_access_iterator_tag,
                                                std::conditional_t<std::ranges::bidirectional_range<Base>,
                                                                   std::bidirectional_iterator_tag,
                                                                   std::forward_iterator_tag>>;
    using value_type = decltype(std::views::take(std::ranges::subrange{current_, end_}, n_));
    using difference_type = std::ranges::range_difference_t<Base>;

    iterator() = default;
    constexpr iterator(iterator<!Const> i)
        requires Const && std::convertible_to<std::ranges::iterator_t<V>, std::ranges::iterator_t<Base>>
                  && std::convertible_to<std::ranges::sentinel_t<V>, std::ranges::sentinel_t<Base>>
        : current_{std::move(i.current_)}, end_{std::move(i.end_)}, n_{i.n_}, missing_{i.missing_}
    {}

    constexpr std::ranges::iterator_t<Base> base() const
    {
        return current_;
    }

    constexpr value_type operator*() const
    {
        assert(current_ != end_);
        return std::views::take(std::ranges::subrange(current_, end_), n_);
    }

    constexpr iterator & operator++()
    {
        assert(current_ != end_);
        missing_ = std::ranges::advance(current_, n_, end_);
        return *this;
    }

    constexpr iterator operator++(int)
    {
        auto tmp = *this;
        ++*this;
        return tmp;
    }

    constexpr iterator & operator--()
        requires std::ranges::bidirectional_range<Base>
    {
        std::ranges::advance(current_, missing_ - n_);
        missing_ = 0;
        return *this;
    }

    constexpr iterator operator--(int)
        requires std::ranges::bidirectional_range<Base>
    {
        auto tmp = *this;
        --*this;
        return tmp;
    }

    constexpr iterator & operator+=(difference_type x)
        requires std::ranges::random_access_range<Base>
    {
        assert(x <= 0 || std::ranges::distance(current_, end_) > n_ * (x - 1));
        if (x > 0)
        {
            std::ranges::advance(current_, n_ * (x - 1));
            missing_ = std::ranges::advance(current_, n_, end_);
        }
        else if (x < 0)
        {
            std::ranges::advance(current_, n_ * x + missing_);
            missing_ = 0;
        }
        return *this;
    }

    constexpr iterator & operator-=(difference_type x)
        requires std::ranges::random_access_range<Base>
    {
        return *this += -x;
    }

    constexpr value_type operator[](difference_type n) const
        requires std::ranges::random_access_range<Base>
    {
        return *(*this + n);
    }

    friend constexpr bool operator==(iterator const & x, iterator const & y)
    {
        return x.current_ == y.current_;
    }

    friend constexpr bool operator==(iterator const & x, std::default_sentinel_t)
    {
        return x.current_ == x.end_;
    }

    friend constexpr bool operator<(iterator const & x, iterator const & y)
        requires std::ranges::random_access_range<Base>
    {
        return x.current_ < y.current_;
    }

    friend constexpr bool operator>(iterator const & x, iterator const & y)
        requires std::ranges::random_access_range<Base>
    {
        return y < x;
    }

    friend constexpr bool operator<=(iterator const & x, iterator const & y)
        requires std::ranges::random_access_range<Base>
    {
        return !(y < x);
    }

    friend constexpr bool operator>=(iterator const & x, iterator const & y)
        requires std::ranges::random_access_range<Base>
    {
        return !(x < y);
    }

#    ifdef __cpp_lib_three_way_comparison
    friend constexpr auto operator<=>(iterator const & x, iterator const & y)
        requires std::ranges::random_access_range<Base> && std::three_way_comparable<std::ranges::iterator_t<Base>>
    {
        return x.current_ <=> y.current_;
    }
#    endif

    friend constexpr iterator operator+(iterator const & i, difference_type n)
        requires std::ranges::random_access_range<Base>
    {
        auto r = i;
        r += n;
        return r;
    }

    friend constexpr iterator operator+(difference_type n, iterator const & i)
        requires std::ranges::random_access_range<Base>
    {
        auto r = i;
        r += n;
        return r;
    }

    friend constexpr iterator operator-(iterator const & i, difference_type n)
        requires std::ranges::random_access_range<Base>
    {
        auto r = i;
        r -= n;
        return r;
    }

    friend constexpr difference_type operator-(iterator const & x, iterator const & y)
        requires std::sized_sentinel_for<std::ranges::sentinel_t<Base>, std::ranges::iterator_t<Base>>
    {
        return (x.current_ - y.current_ + x.missing_ - y.missing_) / x.n_;
    }

    friend constexpr difference_type operator-(std::default_sentinel_t, iterator const & x)
        requires std::sized_sentinel_for<std::ranges::sentinel_t<Base>, std::ranges::iterator_t<Base>>
    {
        return seqan::stl::detail::chunk::div_ceil(x.end_ - x.current_, x.n_);
    }

    friend constexpr difference_type operator-(iterator const & x, std::default_sentinel_t y)
        requires std::sized_sentinel_for<std::ranges::sentinel_t<Base>, std::ranges::iterator_t<Base>>
    {
        return -(y - x);
    }
};

struct chunk_fn
{
    template <typename Difference>
    constexpr auto operator()(Difference n) const
    {
        return seqan::stl::detail::adaptor_from_functor{*this, n};
    }

    template <seqan::stl::ranges::viewable_range Range, typename Difference = std::ranges::range_difference_t<Range>>
    constexpr auto operator()(Range && range, std::type_identity_t<Difference> n) const
    {
        return chunk_view{std::forward<Range>(range), n};
    }
};

} // namespace seqan::stl::ranges

namespace seqan::stl::views
{

inline constexpr auto chunk = seqan::stl::ranges::chunk_fn{};

} // namespace seqan::stl::views

#endif // ifdef __cpp_lib_ranges_chunk

#endif // SEQAN_STD_CHUNK_VIEW
