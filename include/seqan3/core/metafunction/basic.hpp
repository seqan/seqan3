// -----------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2019, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2019, Knut Reinert & MPI für molekulare Genetik
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at: https://github.com/seqan/seqan3/blob/master/LICENSE.md
// -----------------------------------------------------------------------------------------------------

/*!\file
 * \brief Provides various metafunctions on generic types.
 * \author Hannes Hauswedell <hannes.hauswedell AT fu-berlin.de>
 */

#pragma once

#include <tuple>
#include <type_traits>

#include <seqan3/core/metafunction/function.hpp>

namespace seqan3
{

/*!\addtogroup metafunction
 * \{
 */
 // ----------------------------------------------------------------------------
 // delete_const
 // ----------------------------------------------------------------------------

/*!\brief Return the input type with the topmost `const` removed [Type metafunction].
 * \tparam t The type to operate on.
 * \attention In contrast to std::remove_const, this function also removes `const` from pointers and references.
 *
 * | t                | delete_const<t>::type |
 * |------------------|-----------------------|
 * | t                | t                     |
 * | t const          | t                     |
 * | t &              | t &                   |
 * | t const &        | t &                   |
 * | t &&             | t &&                  |
 * | t const &&       | t &&                  |
 * | t *              | t                     |
 * | t const *        | t *                   |
 * | t const * const  | t const *             |
 * | t volatile       | t volatile            |
 * | t const volatile | t volatile            |
 */
template <typename t>
struct delete_const
{
    //!\brief The return type is the input type with the topmost `const` stripped.
    using type = t;
};

//!\cond
template <typename t>
struct delete_const<t const>
{
    using type = t;
};

template <typename t>
struct delete_const<t const *>
{
    using type = t *;
};

template <typename t>
struct delete_const<t const &>
{
    using type = t &;
};

template <typename t>
struct delete_const<t const &&>
{
    using type = t &&;
};
//!\endcond

/*!\brief Return the input type with topmost `const` removed [Type metafunction, shortcut].
 * \tparam t The type to operate on.
 * \attention In contrast to std::remove_const, this function also removes `const` from pointers and references.
 */
template <typename t>
using delete_const_t = typename delete_const<t>::type;

// ----------------------------------------------------------------------------
// remove_cvref_t
// ----------------------------------------------------------------------------

/*!\brief Return the input type with `const`, `volatile` and references removed [Type metafunction].
 * \tparam t The type to operate on.
 */
template <typename t>
using remove_cvref_t = std::remove_cv_t<std::remove_reference_t<t>>;

// ----------------------------------------------------------------------------
// remove_rvalue_reference
// ----------------------------------------------------------------------------

/*!\brief Return the input type with `&&` removed, but lvalue references preserved [Type metafunction].
 * \tparam t The type to operate on.
 * \see seqan3::remove_rvalue_reference_t
 */
template <typename t>
struct remove_rvalue_reference
{
    //!\brief The return type is the input type with any `&&` stripped.
    using type = std::conditional_t<std::is_rvalue_reference_v<t>, std::remove_reference_t<t>, t>;
};

/*!\brief Return the input type with `&&` removed, but lvalue references preserved [Type metafunction, shortcut].
 * \tparam t The type to operate on.
 * \see seqan3::remove_rvalue_reference
 */
template <typename t>
using remove_rvalue_reference_t = typename remove_rvalue_reference<t>::type;

// ----------------------------------------------------------------------------
// is_constexpr_default_constructible
// ----------------------------------------------------------------------------

/*!\brief Whether a type std::is_default_constructible in `constexpr`-context [UnaryTypeTrait].
 * \tparam t The type to operate on.
 */
template <typename t>
struct is_constexpr_default_constructible : std::false_type
{};

/*!\brief Whether a type std::is_default_constructible in `constexpr`-context [UnaryTypeTrait specialisation].
 * \tparam t A type that std::is_default_constructible.
 */
template <typename t>
//!\cond
    requires std::is_default_constructible_v<t>
//!\endcond
struct is_constexpr_default_constructible<t> : std::integral_constant<bool, SEQAN3_IS_CONSTEXPR(t{})>
{};

/*!\brief Whether a type std::is_default_constructible in `constexpr`-context [UnaryTypeTrait shortcut].
 * \tparam t The type to operate on.
 * \see seqan3::is_constexpr_default_constructible
 */
template <typename t>
inline constexpr bool is_constexpr_default_constructible_v = is_constexpr_default_constructible<t>::value;

//!\}
} // namespace seqan3

namespace seqan3::detail
{

/*!\addtogroup metafunction
 * \{
 */

// ----------------------------------------------------------------------------
// deferred_type
// ----------------------------------------------------------------------------

/*!\brief Return the type identity; further arguments are ignored, but can make this type dependent if they are
 *        [TransformationTrait].
 * \tparam t The type to operate on.
 * \tparam dependent_ts Any provided types are ignored.
 */
template <typename t, typename ...dependent_ts>
struct deferred_type
{
    //!\brief The type identity.
    using type = t;
};

/*!\brief Return the type identity; further arguments are ignored, but can make this type dependent if they are
 *        [TransformationTrait shortcut].
 * \tparam t The type to operate on.
 * \tparam dependent_ts Any provided types are ignored.
 * \see seqan3::deferred_type
 */
template <typename t, typename ...dependent_ts>
using deferred_type_t = typename deferred_type<t, dependent_ts...>::type;

// ----------------------------------------------------------------------------
// remove_cvref_t
// ----------------------------------------------------------------------------

//!\brief Return the type of std::ignore with `const`, `volatile` and references removed [Type metafunction].
using ignore_t = remove_cvref_t<decltype(std::ignore)>;

/*!\brief Return whether the input type with `const`, `volatile` and references removed is std::ignore's type.
 * [Value metafunction].
 * \tparam t The type to operate on.
 */
template <typename t>
constexpr bool decays_to_ignore_v = std::is_same_v<remove_cvref_t<t>, ignore_t>;

//!\}

} // namespace seqan3::detail
