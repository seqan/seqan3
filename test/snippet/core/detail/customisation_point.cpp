// SPDX-FileCopyrightText: 2006-2025 Knut Reinert & Freie Universität Berlin
// SPDX-FileCopyrightText: 2016-2025 Knut Reinert & MPI für molekulare Genetik
// SPDX-License-Identifier: CC0-1.0

#include <concepts>
#include <vector>

#include <seqan3/core/detail/customisation_point.hpp>

//! [CPO Definition]
namespace seqan3::detail::adl_only
{
// Poison-pill overload to prevent non-ADL forms of unqualified lookup.
template <typename... args_t>
void begin(args_t...) = delete;

struct begin_cpo : public detail::customisation_point_object<begin_cpo, 1>
{
    using base_t = detail::customisation_point_object<begin_cpo, 1>;
    // Only this class is allowed to import the constructors from base_t. (CRTP safety idiom)
    using base_t::base_t;

    // range.begin(), member access
    //! [SEQAN3_CPO_OVERLOAD]
    template <typename range_t>
        requires true // further constraints
    static constexpr auto SEQAN3_CPO_OVERLOAD(seqan3::detail::priority_tag<1>, range_t && range)(
        /*return*/ std::forward<range_t>(range).begin() /*;*/
    );
    //! [SEQAN3_CPO_OVERLOAD]

    // begin(range), ADL access
    template <typename range_t>
    static constexpr auto SEQAN3_CPO_OVERLOAD(seqan3::detail::priority_tag<0>, range_t && range)(
        /*return*/ begin(std::forward<range_t>(range)) /*;*/
    );
};

} // namespace seqan3::detail::adl_only
//! [CPO Definition]

//! [CPO Instance]
namespace seqan3
{
// CPO is a normal function object that can be called via seqan3::begin(...)
inline constexpr auto begin = detail::adl_only::begin_cpo{};
} // namespace seqan3
//! [CPO Instance]

//! [ADL Definition]
namespace other_library
{
struct foo
{
    friend int begin(foo const &) // ADL begin, as friend
    {
        return 0;
    }
};
} // namespace other_library
//! [ADL Definition]

//! [CPO Member overload]
// seqan3::begin CPO that will call the "begin" member function
std::vector<int> vec{};

static_assert(std::same_as<decltype(seqan3::begin(vec)), decltype(vec.begin())>); // same iterator type

static_assert(noexcept(vec.begin()));                                 // is noexcept
static_assert(noexcept(seqan3::begin(vec)) == noexcept(vec.begin())); // perfect noexcept-forwarding
//! [CPO Member overload]

//! [CPO ADL overload]
// seqan3::begin CPO that will call the "begin" function per ADL
other_library::foo foo{};

static_assert(std::same_as<decltype(seqan3::begin(foo)), decltype(begin(foo))>); // same value type

static_assert(!noexcept(begin(foo)));                                // isn't noexcept
static_assert(noexcept(seqan3::begin(foo)) == noexcept(begin(foo))); // perfect noexcept-forwarding
//! [CPO ADL overload]

//! [CPO SFINAE friendly]
auto cpo_is_sfinae_friendly(...) -> void;

template <typename range_t>
auto cpo_is_sfinae_friendly(range_t && range) -> decltype(seqan3::begin(range));

// seqan3::begin itself is SFINAE friendly, i.e. no-hard compiler errors, if no cpo overload matches
static_assert(std::same_as<decltype(cpo_is_sfinae_friendly(0)), void>);

static_assert(std::same_as<decltype(cpo_is_sfinae_friendly(vec)), decltype(vec.begin())>);
//! [CPO SFINAE friendly]
