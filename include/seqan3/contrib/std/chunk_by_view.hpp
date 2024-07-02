// SPDX-FileCopyrightText: 2006-2024 Knut Reinert & Freie Universität Berlin
// SPDX-FileCopyrightText: 2016-2024 Knut Reinert & MPI für molekulare Genetik
// SPDX-License-Identifier: BSD-3-Clause

/*!\file
 * \author Enrico Seiler <enrico.seiler AT fu-berlin.de>
 * \brief Provides seqan::stl::views::chunk_by.
 */

// File might be included from multiple libraries.
#ifndef SEQAN_STD_CHUNK_BY_VIEW
#define SEQAN_STD_CHUNK_BY_VIEW

#include <algorithm>
#include <cassert>
#include <ranges>

#ifdef __cpp_lib_ranges_chunk_by

namespace seqan::stl::views
{

using std::ranges::views::chunk_by;

} // namespace seqan::stl::views

#else

#    include "all_view.hpp"
#    include "concepts.hpp"
#    include "detail/adaptor_from_functor.hpp"
#    include "detail/compiler_definitions.hpp"
#    include "detail/movable_box.hpp"
#    include "detail/non_propagating_cache.hpp"

namespace seqan::stl::ranges
{

template <std::ranges::forward_range V,
          std::indirect_binary_predicate<std::ranges::iterator_t<V>, std::ranges::iterator_t<V>> Pred>
    requires std::ranges::view<V> && std::is_object_v<Pred>
class chunk_by_view : public std::ranges::view_interface<chunk_by_view<V, Pred>>
{
private:
    V base_ = V();
    seqan::stl::detail::movable_box<Pred> pred_;
    class iterator;

    seqan::stl::detail::non_propagating_cache<std::ranges::iterator_t<V>> current_;

    constexpr std::ranges::iterator_t<V> find_next(std::ranges::iterator_t<V> current)
    {
        assert(pred_.has_value());
        auto pred = [this]<typename lhs_t, typename rhs_t>(lhs_t && lhs, rhs_t && rhs)
        {
            // Call pred with two arguments, return logical not.
            return !bool((*pred_)(std::forward<lhs_t>(lhs), std::forward<rhs_t>(rhs)));
        };
        auto it = std::ranges::adjacent_find(current, std::ranges::end(base_), pred);
        return std::ranges::next(it, 1, std::ranges::end(base_));
    }

    constexpr std::ranges::iterator_t<V> find_prev(std::ranges::iterator_t<V> current)
        requires std::ranges::bidirectional_range<V>
    {
        assert(pred_.has_value());
        auto pred = [this]<typename lhs_t, typename rhs_t>(lhs_t && lhs, rhs_t && rhs)
        {
            // Call pred with two (reversed) arguments, return logical not.
            return !bool((*pred_)(std::forward<rhs_t>(rhs), std::forward<lhs_t>(lhs)));
        };
        auto rbegin = std::make_reverse_iterator(current);
        auto rend = std::make_reverse_iterator(std::ranges::begin(base_));
        assert(rbegin != rend);
        auto it = std::ranges::adjacent_find(rbegin, rend, pred).base();
        return std::ranges::prev(it, 1, std::ranges::begin(base_));
    }

public:
    chunk_by_view()
        requires std::default_initializable<V> && std::default_initializable<Pred>
    = default;

    constexpr explicit chunk_by_view(V base, Pred pred) : base_{std::move(base)}, pred_{std::move(pred)}
    {}

    constexpr V base() const &
        requires std::copy_constructible<V>
    {
        return base_;
    }
    constexpr V base() &&
    {
        return std::move(base_);
    }

    constexpr Pred const & pred() const
    {
        return *pred_;
    }

    constexpr iterator begin()
    {
        assert(pred_.has_value());
        std::ranges::iterator_t<V> it;
        if (current_.has_value())
        {
            it = current_.value();
        }
        else
        {
            it = find_next(std::ranges::begin(base_));
            current_ = it;
        }

        return iterator{*this, std::ranges::begin(base_), it};
    }
    constexpr auto end()
    {
        if constexpr (std::ranges::common_range<V>)
        {
            return iterator{*this, std::ranges::end(base_), std::ranges::end(base_)};
        }
        else
        {
            return std::default_sentinel;
        }
    }
};

template <class R, class Pred>
chunk_by_view(R &&, Pred) -> chunk_by_view<seqan::stl::views::all_t<R>, Pred>;

template <std::ranges::forward_range V,
          std::indirect_binary_predicate<std::ranges::iterator_t<V>, std::ranges::iterator_t<V>> Pred>
    requires std::ranges::view<V> && std::is_object_v<Pred>
class chunk_by_view<V, Pred>::iterator
{
private:
    chunk_by_view * parent_ = nullptr;
    std::ranges::iterator_t<V> current_ = std::ranges::iterator_t<V>{};
    std::ranges::iterator_t<V> next_ = std::ranges::iterator_t<V>{};

    constexpr iterator(chunk_by_view & parent, std::ranges::iterator_t<V> current, std::ranges::iterator_t<V> next) :
        parent_{std::addressof(parent)},
        current_{std::move(current)},
        next_{std::move(next)}
    {}

    friend chunk_by_view;

public:
    using value_type = std::ranges::subrange<std::ranges::iterator_t<V>>;
    using difference_type = std::ranges::range_difference_t<V>;
    using iterator_category = std::input_iterator_tag;
    using iterator_concept = std::conditional_t<std::ranges::bidirectional_range<V>, //
                                                std::bidirectional_iterator_tag,
                                                std::forward_iterator_tag>;

    iterator() = default;

    constexpr value_type operator*() const
    {
        assert(current_ != next_);
        return std::ranges::subrange{current_, next_};
    }
    constexpr iterator & operator++()
    {
        assert(current_ != next_);
        current_ = next_;
        next_ = parent_->find_next(current_);
        return *this;
    }
    constexpr iterator operator++(int)
    {
        auto tmp = *this;
        ++*this;
        return tmp;
    }

    constexpr iterator & operator--()
        requires std::ranges::bidirectional_range<V>
    {
        next_ = current_;
        current_ = parent_->find_prev(next_);
        return *this;
    }
    constexpr iterator operator--(int)
        requires std::ranges::bidirectional_range<V>
    {
        auto tmp = *this;
        --*this;
        return tmp;
    }

    friend constexpr bool operator==(iterator const & x, iterator const & y)
    {
        return x.current_ == y.current_;
    }
    friend constexpr bool operator==(iterator const & x, std::default_sentinel_t)
    {
        return x.current_ == x.next_;
    }
};

struct chunk_by_fn
{
    template <typename Predicate>
    constexpr auto operator()(Predicate && pred) const
    {
        return seqan::stl::detail::adaptor_from_functor{*this, std::forward<Predicate>(pred)};
    }

    template <seqan::stl::ranges::viewable_range Range, typename Predicate>
    constexpr auto operator()(Range && range, Predicate && pred) const
    {
        return chunk_by_view{std::forward<Range>(range), std::forward<Predicate>(pred)};
    }
};

} // namespace seqan::stl::ranges

namespace seqan::stl::views
{

inline constexpr auto chunk_by = seqan::stl::ranges::chunk_by_fn{};

} // namespace seqan::stl::views

#endif // ifdef __cpp_lib_ranges_chunk_by

#endif // SEQAN_STD_CHUNK_BY_VIEW
