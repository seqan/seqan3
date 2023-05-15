// -----------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2023, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2023, Knut Reinert & MPI für molekulare Genetik
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at: https://github.com/seqan/seqan-std/blob/main/LICENSE
// -----------------------------------------------------------------------------------------------------

/*!\file
 * \author Enrico Seiler <enrico.seiler AT fu-berlin.de>
 * \brief Provides seqan::std::views::{all, all_t}, and seqan::std::ranges::owning_view.
 */

// File might be included from multiple libraries.
#ifndef SEQAN_STD_ALL_VIEW
#define SEQAN_STD_ALL_VIEW

#include <ranges>

#if __cpp_lib_ranges >= 202110L

namespace seqan::std::ranges
{

using ::std::ranges::owning_view;

} // namespace seqan::std::ranges

namespace seqan::std::views
{

using ::std::ranges::views::all;
using ::std::ranges::views::all_t;

} // namespace seqan::std::views
#else
#    include "concepts.hpp"
#    include "detail/adaptor_base.hpp"
#    include "detail/exposition_only.hpp"

namespace seqan::std::ranges
{

/*!\brief A move-only view that takes unique ownership of a range.
 * \sa https://en.cppreference.com/w/cpp/ranges/owning_view
 */
template <::std::ranges::range rng_t>
    requires ::std::movable<rng_t> && (!seqan::std::detail::is_initializer_list<::std::remove_cvref_t<rng_t>>)
class owning_view : public ::std::ranges::view_interface<owning_view<rng_t>>
{
private:
    //!\brief The stored range.
    rng_t rng = rng_t();

public:
    owning_view()
        requires ::std::default_initializable<rng_t>
    = default;                                             //!< Defaulted.
    owning_view(owning_view const &) = delete;             //!< Deleted.
    owning_view & operator=(owning_view const &) = delete; //!< Deleted.
    owning_view(owning_view &&) = default;                 //!< Defaulted.
    owning_view & operator=(owning_view &&) = default;     //!< Defaulted.
    ~owning_view() = default;                              //!< Defaulted.

    //!\brief Move construct from a range.
    constexpr owning_view(rng_t && r) : rng(::std::move(r))
    {}

    //!\brief Returns the range.
    constexpr rng_t & base() & noexcept
    {
        return rng;
    }

    //!\overload
    constexpr rng_t const & base() const & noexcept
    {
        return rng;
    }

    //!\overload
    constexpr rng_t && base() && noexcept
    {
        return ::std::move(rng);
    }

    //!\overload
    constexpr rng_t const && base() const && noexcept
    {
        return ::std::move(rng);
    }

    //!\brief Return the begin iterator of the range.
    constexpr ::std::ranges::iterator_t<rng_t> begin()
    {
        return ::std::ranges::begin(rng);
    }

    //!\overload
    constexpr auto begin() const
        requires ::std::ranges::range<rng_t const>
    {
        return ::std::ranges::begin(rng);
    }

    //!\brief Return the end iterator of the range.
    constexpr ::std::ranges::sentinel_t<rng_t> end()
    {
        return ::std::ranges::end(rng);
    }

    //!\overload
    constexpr auto end() const
        requires ::std::ranges::range<rng_t const>
    {
        return ::std::ranges::end(rng);
    }

    //!\brief Checks whether the range is empty.
    constexpr bool empty()
        requires requires { ::std::ranges::empty(rng); }
    {
        return ::std::ranges::empty(rng);
    }

    //!\overload
    constexpr bool empty() const
        requires requires { ::std::ranges::empty(rng); }
    {
        return ::std::ranges::empty(rng);
    }

    //!\brief Returns the size of the range.
    constexpr auto size()
        requires ::std::ranges::sized_range<rng_t>
    {
        return ::std::ranges::size(rng);
    }

    //!\overload
    constexpr auto size() const
        requires ::std::ranges::sized_range<rng_t const>
    {
        return ::std::ranges::size(rng);
    }

    //!\brief Returns the raw data pointer of the range.
    constexpr auto data()
        requires ::std::ranges::contiguous_range<rng_t>
    {
        return ::std::ranges::data(rng);
    }

    //!\overload
    constexpr auto data() const
        requires ::std::ranges::contiguous_range<rng_t const>
    {
        return ::std::ranges::data(rng);
    }
};

/*!\brief The functor for seqan::std::views::all.
 */
class all_fn : public seqan::std::detail::adaptor_base<all_fn>
{
private:
    //!\brief Befriend the base class.
    friend seqan::std::detail::adaptor_base<all_fn>;

    //!\brief Checks whether a type is a view.
    template <typename t>
    static constexpr bool decays_to_view = ::std::ranges::view<::std::decay_t<t>>;

    //!\brief Checks whether a type could be used for ::std::ranges::ref_view.
    template <typename t>
    static constexpr bool valid_for_ref_view = requires { ::std::ranges::ref_view(::std::declval<t>()); };

    //!\brief Checks whether a type could be used for seqan3::detail::owning_view.
    template <typename t>
    static constexpr bool valid_for_owning_view = requires { owning_view(::std::declval<t>()); };

public:
    using seqan::std::detail::adaptor_base<all_fn>::adaptor_base;

    /*!\brief Returns a view that includes all elements of the range argument.
     * \sa https://en.cppreference.com/w/cpp/ranges/all_view
     * \details
     * This implements the new C++20 behaviour that is only available with gcc12 and newer.
     * In contrast to the old ::std::views::all, rvalue ranges can be bound.
     * \returns
     *   * `rng` if `rng` is a view.
     *   * A ::std::ranges::ref_view of `rng` if that expression is valid.
     *   * Otherwise, a seqan3::detail::owning_view of `rng`.
     */
    template <seqan::std::ranges::viewable_range rng_t>
        requires decays_to_view<rng_t> || valid_for_ref_view<rng_t> || valid_for_owning_view<rng_t>
    static auto impl(rng_t && rng)
    {
        if constexpr (decays_to_view<rng_t>)
            return ::std::forward<rng_t>(rng);
        else if constexpr (valid_for_ref_view<rng_t>)
            return ::std::ranges::ref_view{::std::forward<rng_t>(rng)};
        else
            return owning_view{::std::forward<rng_t>(rng)};
    }
};

} // namespace seqan::std::ranges

namespace seqan::std::views
{

/*!\copydoc all_fn::impl
 */
inline constexpr auto all = seqan::std::ranges::all_fn{};

/*!\brief Returns the type that results from appying seqan3::detail::all to a range.
 */
template <seqan::std::ranges::viewable_range rng_t>
using all_t = decltype(all(::std::declval<rng_t>()));

} // namespace seqan::std::views

#endif // __cpp_lib_ranges >= 202110L

#endif // SEQAN_STD_ALL_VIEW
