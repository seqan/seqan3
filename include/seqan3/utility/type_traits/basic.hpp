// SPDX-FileCopyrightText: 2006-2025 Knut Reinert & Freie Universität Berlin
// SPDX-FileCopyrightText: 2016-2025 Knut Reinert & MPI für molekulare Genetik
// SPDX-License-Identifier: BSD-3-Clause

/*!\file
 * \brief Provides various type traits on generic types.
 * \author Hannes Hauswedell <hannes.hauswedell AT fu-berlin.de>
 */

#pragma once

#include <tuple>
#include <type_traits>

#include <seqan3/core/platform.hpp>

// ----------------------------------------------------------------------------
// is_constexpr
// ----------------------------------------------------------------------------

//!\cond DEV
/*!\brief Returns true if the expression passed to this macro can be evaluated at compile time, false otherwise.
 * \ingroup utility_type_traits
 * \returns true or false.
 */
#if defined(__clang__) // https://github.com/llvm/llvm-project/issues/58078
#    define SEQAN3_IS_CONSTEXPR(...) true
#else
#    define SEQAN3_IS_CONSTEXPR(...) std::integral_constant<bool, __builtin_constant_p(((void)__VA_ARGS__, 0))>::value
#endif
//!\endcond

namespace seqan3
{

// ----------------------------------------------------------------------------
// remove_rvalue_reference
// ----------------------------------------------------------------------------

/*!\brief Return the input type with `&&` removed, but lvalue references preserved.
 * \ingroup utility_type_traits
 * \implements seqan3::transformation_trait
 * \tparam t The type to operate on.
 * \see seqan3::remove_rvalue_reference_t
 */
template <typename t>
struct remove_rvalue_reference
{
    //!\brief The return type is the input type with any `&&` stripped.
    using type = std::conditional_t<std::is_rvalue_reference_v<t>, std::remove_reference_t<t>, t>;
};

/*!\brief Return the input type with `&&` removed, but lvalue references preserved (transformation_trait shortcut).
 * \tparam t The type to operate on.
 * \see seqan3::remove_rvalue_reference
 * \relates seqan3::remove_rvalue_reference
 */
template <typename t>
using remove_rvalue_reference_t = typename remove_rvalue_reference<t>::type;

// ----------------------------------------------------------------------------
// is_constexpr_default_constructible
// ----------------------------------------------------------------------------

/*!\brief Whether a type std::is_default_constructible in `constexpr`-context (unary_type_trait specialisation).
 * \ingroup utility_type_traits
 * \implements seqan3::unary_type_trait
 * \tparam t The type to operate on.
 */
template <typename t>
struct is_constexpr_default_constructible : std::false_type
{};

//!\cond
template <typename t>
    requires std::is_default_constructible_v<t>
struct is_constexpr_default_constructible<t> : std::integral_constant<bool, SEQAN3_IS_CONSTEXPR(t{})>
{};
//!\endcond

/*!\brief Whether a type std::is_default_constructible in `constexpr`-context (unary_type_trait shortcut).
 * \tparam t The type to operate on.
 * \relates seqan3::is_constexpr_default_constructible
 */
template <typename t>
inline constexpr bool is_constexpr_default_constructible_v = is_constexpr_default_constructible<t>::value;

} // namespace seqan3

namespace seqan3::detail
{

// ----------------------------------------------------------------------------
// deferred_type
// ----------------------------------------------------------------------------

/*!\brief Return the type identity; further arguments are ignored, but can make this type dependent if they are.
 * \ingroup utility_type_traits
 * \implements seqan3::transformation_trait
 * \tparam t The type to operate on.
 * \tparam dependent_ts Any provided types are ignored.
 * \see seqan3::detail::deferred_type_t
 */
template <typename t, typename... dependent_ts>
struct deferred_type
{
    //!\brief The type identity.
    using type = t;
};

/*!\brief Return the type identity; further arguments are ignored, but can make this type dependent if they are
 *        (transformation_trait shortcut).
 * \ingroup utility_type_traits
 * \tparam t The type to operate on.
 * \tparam dependent_ts Any provided types are ignored.
 * \see seqan3::detail::deferred_type
 */
template <typename t, typename... dependent_ts>
using deferred_type_t = typename deferred_type<t, dependent_ts...>::type;

// ----------------------------------------------------------------------------
// ignore_t
// ----------------------------------------------------------------------------

//!\brief Return the type of std::ignore with `const`, `volatile` and references removed (type trait).
//!\ingroup utility_type_traits
using ignore_t = std::remove_cvref_t<decltype(std::ignore)>;

/*!\brief Return whether the input type with `const`, `volatile` and references removed is std::ignore's type.
 *        (type trait).
 * \ingroup utility_type_traits
 * \tparam t The type to operate on.
 */
template <typename t>
constexpr bool decays_to_ignore_v = std::is_same_v<std::remove_cvref_t<t>, ignore_t>;

} // namespace seqan3::detail

// ----------------------------------------------------------------------------
// SEQAN3_IS_SAME
// ----------------------------------------------------------------------------

//!\cond DEV
/*!\brief A macro that behaves like std::is_same_v, except that it doesn't need to instantiate the template on GCC and
 *        Clang.
 * \ingroup utility_type_traits
 */
#if defined(__clang__)
#    define SEQAN3_IS_SAME(...) __is_same(__VA_ARGS__)
#elif defined(__GNUC__)
#    define SEQAN3_IS_SAME(...) __is_same_as(__VA_ARGS__)
#else
#    define SEQAN3_IS_SAME(...) std::is_same_v<__VA_ARGS__>
#endif
//!\endcond
