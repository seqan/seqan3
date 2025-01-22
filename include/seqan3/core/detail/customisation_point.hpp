// SPDX-FileCopyrightText: 2006-2025 Knut Reinert & Freie Universität Berlin
// SPDX-FileCopyrightText: 2016-2025 Knut Reinert & MPI für molekulare Genetik
// SPDX-License-Identifier: BSD-3-Clause

/*!\file
 * \author Hannes Hauswedell <hannes.hauswedell AT fu-berlin.de>
 * \brief Helper utilities for defining customisation point objects (CPOs).
 */

#pragma once

#include <seqan3/core/platform.hpp>

namespace seqan3::detail
{

// ============================================================================
// priority_tag
// ============================================================================

//!\brief A tag that allows controlled overload resolution via implicit base conversion rules.
//!\ingroup core
template <size_t I>
struct priority_tag
    //!\cond
    // Doxygen fail
    : priority_tag<I - 1>
//!\endcond
{};

//!\brief Recursion anchor for seqan3::detail::priority_tag.
template <>
struct priority_tag<0>
{};

} // namespace seqan3::detail

// ============================================================================
// SEQAN3_CPO_OVERLOAD
// ============================================================================
#if SEQAN3_DOXYGEN_ONLY(1) 0
/*!\brief A macro helper for #SEQAN3_CPO_OVERLOAD.
 * \cond DEV
 * \ingroup core
 * \endcond DEV
 * \details
 * Please note that in order to allow a semicolon at the end when using this macro, i.e.
 *
 * ```cpp
 *  SEQAN3_CPO_OVERLOAD_BODY
 *  (
 *      ...
 *  );
 * ```
 *
 * we need some expression at the end of the macro that can have a semicolon appended.
 * We chose `static_assert(true)` for that reason.
 */
#    define SEQAN3_CPO_OVERLOAD_BODY(...)                                                                              \
        noexcept(auto)                                                                                                 \
        {                                                                                                              \
            return __VA_ARGS__;                                                                                        \
        }
#else // ^^^ (simplified) doxygen version / normal definition vvv
#    define SEQAN3_CPO_OVERLOAD_BODY(...)                                                                              \
        noexcept(noexcept(__VA_ARGS__))->decltype(__VA_ARGS__)                                                         \
        {                                                                                                              \
            return __VA_ARGS__;                                                                                        \
        }                                                                                                              \
        static_assert(true)
#endif

/*!\brief A macro that helps to define a seqan3::detail::customisation_point_object.
 * \cond DEV
 * \ingroup core
 * \endcond DEV
 * \details
 * Expands to a function definition with the name `cpo_overload`.
 *
 * It puts the given expression via #SEQAN3_CPO_OVERLOAD_BODY as a single return statement in the function body, the
 * noexcept declaration and requires declaration.
 *
 * \snippet test/snippet/core/detail/customisation_point.cpp SEQAN3_CPO_OVERLOAD
 *
 * expands to something similar to
 *
 * ```cpp
 * template <typename range_t>
 *     requires true // further constraints
 * static constexpr auto cpo_overload(seqan3::detail::priority_tag<1>, range_t && range)
 *
 *     noexcept(noexcept(std::forward<range_t>(range).begin()))
 *
 *     requires requires()
 *     {
 *         {std::forward<range_t>(range).begin()};
 *     }
 *
 * {
 *     return std::forward<range_t>(range).begin();
 * }
 * ```
 */
#define SEQAN3_CPO_OVERLOAD(...) cpo_overload(__VA_ARGS__) SEQAN3_CPO_OVERLOAD_BODY

namespace seqan3::detail
{
/*!\brief A CRTP base-class that defines a customisation_point_object (CPO).
 * \ingroup core
 * \tparam derived_t CRTP derived type.
 * \tparam max_priority The (highest) start priority the CPO overloads will be checked from.
 *
 * You can define a CPO similar to `decltype(std::ranges::begin)` that accepts `range.begin()` (Member) or
 * `begin(range)` (ADL):
 *
 * \snippet test/snippet/core/detail/customisation_point.cpp CPO Definition
 *
 * and instantiate it like std::ranges::begin:
 *
 * \snippet test/snippet/core/detail/customisation_point.cpp CPO Instance
 *
 * Any type that has a member function with the name "begin" will now work (a real definition of std::ranges::begin
 * would, of course, further constraint the overload):
 *
 * \snippet test/snippet/core/detail/customisation_point.cpp CPO Member overload
 *
 * Furthermore, any unqualified function call "begin":
 *
 * \snippet test/snippet/core/detail/customisation_point.cpp CPO ADL overload
 *
 * seqan3::detail::customisation_point_object is designed to be SFINAE friendly:
 *
 * \snippet test/snippet/core/detail/customisation_point.cpp CPO SFINAE friendly
 *
 */
template <typename derived_t, unsigned max_priority>
struct customisation_point_object
{
private:
    //!\brief Allow `derived_t` to inherit the constructors of this CRTP-class.
    friend derived_t;

    constexpr customisation_point_object() = default;                                               //!< Defaulted.
    constexpr customisation_point_object(customisation_point_object &&) = default;                  //!< Defaulted.
    constexpr customisation_point_object(customisation_point_object const &) = default;             //!< Defaulted.
    constexpr customisation_point_object & operator=(customisation_point_object &&) = default;      //!< Defaulted.
    constexpr customisation_point_object & operator=(customisation_point_object const &) = default; //!< Defaulted.

public:
    /*!\brief SFINAE-friendly call-operator of this seqan3::detail::customisation_point_object instance.
     * \details
     * This operator implements the actual CPO overload resolution. Overload resolution will try each base class of
     * seqan3::detail::priority_tag\<max_priority>, seqan3::detail::priority_tag\<max_priority - 1>, ...,
     * seqan3::detail::priority_tag<0> as first argument to `derived_t::cpo_overload`. That means a high priority in
     * seqan3::detail::priority_tag will be evaluated first and one can define an order which overload should be
     * prioritised if multiple overloads match.
     *
     * It perfectly forwards the result and noexcept-property of the `derived_t::cpo_overload`.
     */
    template <typename... args_t, typename derived_type = derived_t /*circumvent incomplete types*/>
    constexpr auto operator()(args_t &&... args) const SEQAN3_CPO_OVERLOAD_BODY(
        /*return*/ derived_type::cpo_overload(priority_tag<max_priority>{}, std::forward<args_t>(args)...) /*;*/
    );
};
} // namespace seqan3::detail
