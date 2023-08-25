// -----------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2023, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2023, Knut Reinert & MPI für molekulare Genetik
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at: https://github.com/seqan/seqan-std/blob/main/LICENSE
// -----------------------------------------------------------------------------------------------------

/*!\file
 * \author Enrico Seiler <enrico.seiler AT fu-berlin.de>
 * \brief Provides seqan::stl::views::join_with.
 */

// File might be included from multiple libraries.
#ifndef SEQAN_STD_JOIN_WITH_VIEW
#define SEQAN_STD_JOIN_WITH_VIEW

#include <ranges>

#ifdef __cpp_lib_ranges_join_with

namespace seqan::stl::views
{

using std::ranges::views::join_with;

} // namespace seqan::stl::views

#else

#    include <variant>

#    include "all_view.hpp"
#    include "concepts.hpp"
#    include "detail/adaptor_from_functor.hpp"
#    include "detail/compiler_definitions.hpp"
#    include "detail/exposition_only.hpp"
#    include "detail/non_propagating_cache.hpp"

namespace seqan::stl::detail::join_with
{

template <class R, class P>
concept compatible_joinable_ranges =
    std::common_with<std::ranges::range_value_t<R>, std::ranges::range_value_t<P>>
    && std::common_reference_with<std::ranges::range_reference_t<R>, std::ranges::range_reference_t<P>>
    && std::common_reference_with<std::ranges::range_rvalue_reference_t<R>, std::ranges::range_rvalue_reference_t<P>>;

template <class R>
concept bidirectional_common = std::ranges::bidirectional_range<R> && std::ranges::common_range<R>;

// A helper to decide whether the join_with view should have a non_propagating_cache member.
template <typename InnerRng>
struct cache_helper
{};

// If the inner range is a reference, we do not need a cache.
template <typename InnerRng>
    requires std::is_reference_v<InnerRng>
struct cache_helper<InnerRng>
{};

// If the inner range is not a reference, we need a cache.
template <typename InnerRng>
    requires (!std::is_reference_v<InnerRng>)
struct cache_helper<InnerRng>
{
    non_propagating_cache<std::remove_cv_t<InnerRng>> inner_;
};

} // namespace seqan::stl::detail::join_with

namespace seqan::stl::ranges
{

template <std::ranges::input_range V, std::ranges::forward_range Pattern>
    requires std::ranges::view<V> && std::ranges::input_range<std::ranges::range_reference_t<V>>
              && std::ranges::view<Pattern>
              && seqan::stl::detail::join_with::compatible_joinable_ranges<std::ranges::range_reference_t<V>, Pattern>
class join_with_view :
    public std::ranges::view_interface<join_with_view<V, Pattern>>,
    private seqan::stl::detail::join_with::cache_helper<std::ranges::range_reference_t<V>>
{
private:
    using InnerRng = std::ranges::range_reference_t<V>;

    V base_{};

    Pattern pattern_{};

    template <bool Const>
    struct iterator;

    template <bool Const>
    struct sentinel;

public:
    join_with_view()
        requires std::default_initializable<V> && std::default_initializable<Pattern>
    = default;

    constexpr join_with_view(V base, Pattern pattern) : base_(std::move(base)), pattern_(std::move(pattern))
    {}

    template <std::ranges::input_range R>
        requires std::constructible_from<V, std::ranges::views::all_t<R>>
                  && std::constructible_from<Pattern, std::ranges::single_view<std::ranges::range_value_t<InnerRng>>>
    constexpr join_with_view(R && r, std::ranges::range_value_t<InnerRng> e) :
        base_(std::ranges::views::all(std::forward<R>(r))),
        pattern_(std::ranges::views::single(std::move(e)))
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

    constexpr auto begin()
    {
        // In essence, simple_view means that the iterators of the non-const view and const range are the same and
        // the sentinels of the non-const view and const range are the same.
        // If std::is_reference_v<InnerRng>, we do not have a cache to store the value. Hence, we do not change
        // any members and can use the const iterator.
        constexpr bool use_const = seqan::stl::detail::simple_view<V> && std::is_reference_v<InnerRng>
                                && seqan::stl::detail::simple_view<Pattern>;
        return iterator<use_const>{*this, std::ranges::begin(base_)};
    }

    constexpr auto begin() const
        requires std::ranges::input_range<V const> && std::ranges::forward_range<Pattern const>
              && std::is_reference_v<std::ranges::range_reference_t<V const>>
    {
        return iterator<true>{*this, std::ranges::begin(base_)};
    }

    constexpr auto end()
    {
        constexpr bool is_simple = seqan::stl::detail::simple_view<V> && seqan::stl::detail::simple_view<Pattern>;
        if constexpr (std::ranges::forward_range<V> && std::is_reference_v<InnerRng>
                      && std::ranges::forward_range<InnerRng> && std::ranges::common_range<V>
                      && std::ranges::common_range<InnerRng>)
            return iterator<is_simple>{*this, std::ranges::end(base_)};
        else
            return sentinel<is_simple>{*this};
    }

    constexpr auto end() const
        requires std::ranges::input_range<V const> && std::ranges::forward_range<Pattern const>
              && std::is_reference_v<std::ranges::range_reference_t<V const>>
    {
        using InnerConstRng = std::ranges::range_reference_t<V const>;
        if constexpr (std::ranges::forward_range<V const> && std::ranges::forward_range<InnerConstRng>
                      && std::ranges::common_range<V const> && std::ranges::common_range<InnerConstRng>)
            return iterator<true>{*this, std::ranges::end(base_)};
        else
            return sentinel<true>{*this};
    }
};

template <class R, class P>
join_with_view(R &&, P &&) -> join_with_view<seqan::stl::views::all_t<R>, seqan::stl::views::all_t<P>>;

template <std::ranges::input_range R>
join_with_view(R &&, std::ranges::range_value_t<std::ranges::range_reference_t<R>>)
    -> join_with_view<seqan::stl::views::all_t<R>,
                      std::ranges::single_view<std::ranges::range_value_t<std::ranges::range_reference_t<R>>>>;

} // namespace seqan::stl::ranges

namespace seqan::stl::detail::join_with
{

template <bool Const, typename V, typename Pattern>
struct helper
{
    using Parent = seqan::stl::detail::maybe_const<Const, seqan::stl::ranges::join_with_view<V, Pattern>>;
    using Base = seqan::stl::detail::maybe_const<Const, V>;
    using InnerBase = std::ranges::range_reference_t<Base>;
    using PatternBase = seqan::stl::detail::maybe_const<Const, Pattern>;

    using OuterIter = std::ranges::iterator_t<Base>;
    using InnerIter = std::ranges::iterator_t<InnerBase>;
    using PatternIter = std::ranges::iterator_t<PatternBase>;

    static constexpr bool ref_is_glvalue = std::is_reference_v<InnerBase>;
};

template <bool Const, typename V, typename Pattern, bool ref_is_glvalue>
struct iterator_category_t;

// If std::is_reference_v<InnerBase> is true, define iterator_category.
template <bool Const, typename V, typename Pattern>
struct iterator_category_t<Const, V, Pattern, true>
{
    using Parent = helper<Const, V, Pattern>::Parent;
    using Base = helper<Const, V, Pattern>::Base;
    using InnerBase = helper<Const, V, Pattern>::InnerBase;
    using PatternBase = helper<Const, V, Pattern>::PatternBase;

    using OuterIter = helper<Const, V, Pattern>::OuterIter;
    using InnerIter = helper<Const, V, Pattern>::InnerIter;
    using PatternIter = helper<Const, V, Pattern>::PatternIter;

    using iterator_category = std::conditional_t<
        !std::is_lvalue_reference_v<
            std::common_reference_t<std::iter_reference_t<InnerIter>, std::iter_reference_t<PatternIter>>>,
        std::input_iterator_tag,
        std::conditional_t<
            std::derived_from<std::bidirectional_iterator_tag,
                              typename std::iterator_traits<OuterIter>::iterator_category>
                && std::derived_from<std::bidirectional_iterator_tag,
                                     typename std::iterator_traits<InnerIter>::iterator_category>
                && std::derived_from<std::bidirectional_iterator_tag,
                                     typename std::iterator_traits<PatternIter>::iterator_category>
                && std::ranges::common_range<InnerBase> && std::ranges::common_range<PatternBase>,
            std::bidirectional_iterator_tag,
            std::conditional_t<std::derived_from<std::forward_iterator_tag,
                                                 typename std::iterator_traits<OuterIter>::iterator_category>
                                   && std::derived_from<std::forward_iterator_tag,
                                                        typename std::iterator_traits<InnerIter>::iterator_category>
                                   && std::derived_from<std::forward_iterator_tag,
                                                        typename std::iterator_traits<PatternIter>::iterator_category>,
                               std::forward_iterator_tag,
                               std::input_iterator_tag>>>;
};

// If std::is_reference_v<InnerBase> is false, there is no iterator_category.
template <bool Const, typename V, typename Pattern>
struct iterator_category_t<Const, V, Pattern, false>
{};

} // namespace seqan::stl::detail::join_with

namespace seqan::stl::ranges
{

template <std::ranges::input_range V, std::ranges::forward_range Pattern>
    requires std::ranges::view<V> && std::ranges::input_range<std::ranges::range_reference_t<V>>
          && std::ranges::view<Pattern>
          && seqan::stl::detail::join_with::compatible_joinable_ranges<std::ranges::range_reference_t<V>, Pattern>
template <bool Const>
class join_with_view<V, Pattern>::iterator :
    public seqan::stl::detail::join_with::
        iterator_category_t<Const, V, Pattern, seqan::stl::detail::join_with::helper<Const, V, Pattern>::ref_is_glvalue>
{
private:
    using helper_t = seqan::stl::detail::join_with::helper<Const, V, Pattern>;
    using Parent = helper_t::Parent;
    using Base = helper_t::Base;
    using InnerBase = helper_t::InnerBase;
    using PatternBase = helper_t::PatternBase;

    using OuterIter = helper_t::OuterIter;
    using InnerIter = helper_t::InnerIter;
    using PatternIter = helper_t::PatternIter;

    static constexpr bool ref_is_glvalue = helper_t::ref_is_glvalue;

    friend class join_with_view<V, Pattern>;

    Parent * parent_{nullptr};

    // clang-format off
SEQAN_STD_NESTED_VISIBILITY
    // clang-format on
    OuterIter outer_it_{};

private:
    std::variant<PatternIter, InnerIter> inner_it_{};

    constexpr iterator(Parent & parent, std::ranges::iterator_t<Base> outer) :
        parent_(std::addressof(parent)),
        outer_it_(std::move(outer))
    {
        if (outer_it_ != std::ranges::end(parent_->base_))
        {
            auto && inner = update_inner(outer_it_);
            inner_it_.template emplace<1>(std::ranges::begin(inner));
            satisfy();
        }
    }

    constexpr auto && update_inner(OuterIter const & x)
    {
        if constexpr (ref_is_glvalue)
            return *x;
        else
            return parent_->inner_.emplace_deref(x);
    }

    constexpr auto && get_inner(OuterIter const & x)
    {
        if constexpr (ref_is_glvalue)
            return *x;
        else
            return *parent_->inner_;
    }

    // Skips over empty inner ranges.
    // https://eel.is/c++draft/range.join.with#iterator-7
    constexpr void satisfy()
    {
        while (true)
        {
            if (inner_it_.index() == 0)
            {
                if (std::get<0>(inner_it_) != std::ranges::end(parent_->pattern_))
                    break;
                auto && inner = update_inner(outer_it_);
                inner_it_.template emplace<1>(std::ranges::begin(inner));
            }
            else
            {
                auto && inner = get_inner(outer_it_);
                if (std::get<1>(inner_it_) != std::ranges::end(inner))
                    break;
                if (++outer_it_ == std::ranges::end(parent_->base_))
                {
                    if constexpr (ref_is_glvalue)
                        inner_it_.template emplace<0>();
                    break;
                }
                inner_it_.template emplace<0>(std::ranges::begin(parent_->pattern_));
            }
        }
    }

public:
    using iterator_concept = std::conditional_t<
        !ref_is_glvalue,
        std::input_iterator_tag,
        std::conditional_t<std::ranges::bidirectional_range<Base>
                               && seqan::stl::detail::join_with::bidirectional_common<InnerBase>
                               && seqan::stl::detail::join_with::bidirectional_common<PatternBase>,
                           std::bidirectional_iterator_tag,
                           std::conditional_t<std::ranges::forward_range<Base> && std::ranges::forward_range<InnerBase>,
                                              std::forward_iterator_tag,
                                              std::input_iterator_tag>>>;

    using value_type = std::common_type_t<std::iter_value_t<InnerIter>, std::iter_value_t<PatternIter>>;
    using difference_type = std::common_type_t<std::iter_difference_t<OuterIter>,
                                               std::iter_difference_t<InnerIter>,
                                               std::iter_difference_t<PatternIter>>;

    iterator()
        requires std::default_initializable<OuterIter>
    = default;

    constexpr iterator(iterator<!Const> i)
        requires Const && std::convertible_to<std::ranges::iterator_t<V>, OuterIter>
                  && std::convertible_to<std::ranges::iterator_t<InnerRng>, InnerIter>
                  && std::convertible_to<std::ranges::iterator_t<Pattern>, PatternIter>
        : outer_it_(std::move(i.outer_it_)), parent_(i.parent_)
    {
        if (i.inner_it_.index() == 0)
            inner_it_.template emplace<0>(std::get<0>(std::move(i.inner_it_)));
        else
            inner_it_.template emplace<1>(std::get<1>(std::move(i.inner_it_)));
    }

    constexpr decltype(auto) operator*() const
    {
        using reference = std::common_reference_t<std::iter_reference_t<InnerIter>, std::iter_reference_t<PatternIter>>;
        return std::visit(
            [](auto & it) -> reference
            {
                return *it;
            },
            inner_it_);
    }

    constexpr iterator & operator++()
    {
        std::visit(
            [](auto & it)
            {
                ++it;
            },
            inner_it_);
        satisfy();
        return *this;
    }

    constexpr void operator++(int)
    {
        ++*this;
    }

    constexpr iterator operator++(int)
        requires ref_is_glvalue && std::forward_iterator<OuterIter> && std::forward_iterator<InnerIter>
    {
        iterator tmp = *this;
        ++*this;
        return tmp;
    }

    constexpr iterator & operator--()
        requires ref_is_glvalue && std::ranges::bidirectional_range<Base>
              && seqan::stl::detail::join_with::bidirectional_common<InnerBase>
              && seqan::stl::detail::join_with::bidirectional_common<PatternBase>
    {
        if (outer_it_ == std::ranges::end(parent_->base_))
        {
            auto && inner = *--outer_it_;
            inner_it_.template emplace<1>(std::ranges::end(inner));
        }

        // Similar to satifsy(). Skips over empty inner ranges when going backwards.
        while (true)
        {
            if (inner_it_.index() == 0)
            {
                auto & it = std::get<0>(inner_it_);
                if (it == std::ranges::begin(parent_->pattern_))
                {
                    auto && inner = *--outer_it_;
                    inner_it_.template emplace<1>(std::ranges::end(inner));
                }
                else
                {
                    break;
                }
            }
            else
            {
                auto & it = std::get<1>(inner_it_);
                auto && inner = *outer_it_;
                if (it == std::ranges::begin(inner))
                {
                    inner_it_.template emplace<0>(std::ranges::end(parent_->pattern_));
                }
                else
                {
                    break;
                }
            }
        }

        std::visit(
            [](auto & it)
            {
                --it;
            },
            inner_it_);
        return *this;
    }

    constexpr iterator operator--(int)
        requires ref_is_glvalue && std::ranges::bidirectional_range<Base>
              && seqan::stl::detail::join_with::bidirectional_common<InnerBase>
              && seqan::stl::detail::join_with::bidirectional_common<PatternBase>
    {
        iterator tmp = *this;
        --*this;
        return tmp;
    }

    friend constexpr bool operator==(iterator const & x, iterator const & y)
        requires ref_is_glvalue && std::equality_comparable<OuterIter> && std::equality_comparable<InnerIter>
    {
        return x.outer_it_ == y.outer_it_ && x.inner_it_ == y.inner_it_;
    }

    friend constexpr decltype(auto) iter_move(iterator const & x)
    {
        using rvalue_reference =
            std::common_reference_t<std::iter_rvalue_reference_t<InnerIter>, std::iter_rvalue_reference_t<PatternIter>>;
        return std::visit<rvalue_reference>(std::ranges::iter_move, x.inner_it_);
    }

    friend constexpr void iter_swap(iterator const & x, iterator const & y)
        requires std::indirectly_swappable<InnerIter, PatternIter>
    {
        std::visit(std::ranges::iter_swap, x.inner_it_, y.inner_it_);
    }
};

template <std::ranges::input_range V, std::ranges::forward_range Pattern>
    requires std::ranges::view<V> && std::ranges::input_range<std::ranges::range_reference_t<V>>
          && std::ranges::view<Pattern>
          && seqan::stl::detail::join_with::compatible_joinable_ranges<std::ranges::range_reference_t<V>, Pattern>
template <bool Const>
class join_with_view<V, Pattern>::sentinel
{
private:
    using Parent = seqan::stl::detail::maybe_const<Const, join_with_view>;
    using Base = seqan::stl::detail::maybe_const<Const, V>;
    std::ranges::sentinel_t<Base> end_ = std::ranges::sentinel_t<Base>();

    friend class join_with_view<V, Pattern>;

    constexpr explicit sentinel(Parent & parent) : end_(std::ranges::end(parent.base_))
    {}

public:
    sentinel() = default;
    constexpr sentinel(sentinel<!Const> s)
        requires Const && std::convertible_to<std::ranges::sentinel_t<V>, std::ranges::sentinel_t<Base>>
        : end_(std::move(s.end_))
    {}

    template <bool OtherConst>
        requires std::sentinel_for<std::ranges::sentinel_t<Base>,
                                   std::ranges::iterator_t<seqan::stl::detail::maybe_const<OtherConst, V>>>
    friend constexpr bool operator==(iterator<OtherConst> const & x, sentinel const & y)
    {
        return x.outer_it_ == y.end_;
    }
};

struct join_with_fn
{
    template <typename Pattern>
    constexpr auto operator()(Pattern && pattern) const
    {
        return seqan::stl::detail::adaptor_from_functor{*this, std::forward<Pattern>(pattern)};
    }

    template <seqan::stl::ranges::viewable_range urng_t, typename Pattern>
    constexpr auto operator()(urng_t && urange, Pattern && pattern) const
    {
        return join_with_view{std::forward<urng_t>(urange), std::forward<Pattern>(pattern)};
    }
};

} // namespace seqan::stl::ranges

namespace seqan::stl::views
{

inline constexpr auto join_with = seqan::stl::ranges::join_with_fn{};

} // namespace seqan::stl::views

#endif // ifdef __cpp_lib_ranges_join_with

#endif // SEQAN_STD_JOIN_WITH_VIEW
