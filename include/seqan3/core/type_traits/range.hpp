// -----------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2020, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2020, Knut Reinert & MPI für molekulare Genetik
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at: https://github.com/seqan/seqan3/blob/master/LICENSE.md
// -----------------------------------------------------------------------------------------------------

/*!\file
 * \brief Provides various transformation traits used by the range module.
 * \author Hannes Hauswedell <hannes.hauswedell AT fu-berlin.de>
 */

#pragma once

#include <type_traits>

#include <seqan3/core/platform.hpp>
#include <seqan3/core/type_traits/pre.hpp>
#include <seqan3/core/type_traits/basic.hpp>
#include <seqan3/core/type_traits/iterator.hpp>
#include <seqan3/std/ranges>
#include <seqan3/std/iterator>

// TODO(h-2): add innermost_reference instead of or addition to innermost_value_type?

//NOTE(h-2): for the range overloads we explicitly forbid that the type is iteratoer
// because some types are actually both (e.g. std::directory_iterator)

namespace seqan3::detail
{

//!\cond
template <typename t>
SEQAN3_CONCEPT has_value_type = requires { typename value_type_t<remove_cvref_t<t>>; };
//!\endcond

} // namespace seqan3::detail

namespace seqan3
{

/*!\addtogroup type_traits
 * \{
 */

// ----------------------------------------------------------------------------
// value_type
// ----------------------------------------------------------------------------

/*!\brief Exposes the `value_type` of another type.
 * \implements seqan3::transformation_trait
 * \tparam t The type you wish to query; must model std::ranges::input_range.
 */
template <std::ranges::input_range rng_t>
//!\cond
    requires (!std::input_or_output_iterator<rng_t>)
//!\endcond
struct value_type<rng_t>
{
    //!\brief Return the value_type member definition from the queried type's iterator.
    using type = value_type_t<std::ranges::iterator_t<rng_t>>;
};

// ----------------------------------------------------------------------------
// reference
// ----------------------------------------------------------------------------

/*!\brief Exposes the `reference` of another type.
 * \implements seqan3::transformation_trait
 * \tparam t The type you wish to query; must model std::ranges::input_range.
 */
template <std::ranges::input_range rng_t>
//!\cond
    requires (!std::input_or_output_iterator<rng_t>)
//!\endcond
struct reference<rng_t>
{
    //!\brief Return the reference member definition from the queried type's iterator.
    using type = reference_t<std::ranges::iterator_t<rng_t>>;
};

// ----------------------------------------------------------------------------
// rvalue_reference
// ----------------------------------------------------------------------------

/*!\brief Exposes the `rvalue_reference` of another type.
 * \implements seqan3::transformation_trait
 * \tparam t The type you wish to query; must model std::ranges::input_range.
 */
template <std::ranges::input_range rng_t>
//!\cond
    requires (!std::input_or_output_iterator<rng_t>)
//!\endcond
struct rvalue_reference<rng_t>
{
    //!\brief Return the rvalue_reference member definition from the queried type's iterator.
    using type = rvalue_reference_t<std::ranges::iterator_t<rng_t>>;
};

// ----------------------------------------------------------------------------
// const_reference
// ----------------------------------------------------------------------------

/*!\brief Exposes the `const_reference` of another type.
 * \implements seqan3::transformation_trait
 * \tparam t The type you wish to query; must model std::ranges::input_range.
 */
template <std::ranges::input_range rng_t>
//!\cond
    requires (!std::input_or_output_iterator<rng_t>)
//!\endcond
struct const_reference<rng_t>
{
    //!\brief Resolves to the reference type of the `const_iterator` of t (not the `const iterator`!).
    using type = std::iter_reference_t<std::ranges::iterator_t<rng_t const>>;
};

// ----------------------------------------------------------------------------
// difference_type
// ----------------------------------------------------------------------------

/*!\brief Exposes the `difference_type` of another type.
 * \implements seqan3::transformation_trait
 * \tparam t The type you wish to query; must model std::ranges::input_range.
 */
template <std::ranges::range rng_t>
//!\cond
    requires (!std::input_or_output_iterator<rng_t>)
//!\endcond
struct difference_type<rng_t>
{
    //!\brief Return the difference_type member definition from the queried type's iterator.
    using type = difference_type_t<std::ranges::iterator_t<rng_t>>;
};

// ----------------------------------------------------------------------------
// size_type
// ----------------------------------------------------------------------------

/*!\brief Exposes the `size_type` of another type.
 * \implements seqan3::transformation_trait
 * \tparam t The type you wish to query; must model std::ranges::sized_range.
 */
template <std::ranges::sized_range rng_t>
//!\cond
    requires (!std::input_or_output_iterator<rng_t>)
//!\endcond
struct size_type<rng_t>
{
    //!\brief Return the size_type as returned by the size function.
    using type = decltype(std::ranges::size(std::declval<rng_t &>()));
};

// ----------------------------------------------------------------------------
// innermost_value_type
// ----------------------------------------------------------------------------

//NOTE(h-2): this could be moved to a separate file, because it also applies to iterators

/*!\brief Recursively determines the `value_type` on containers and/or iterators.
 * \implements seqan3::transformation_trait
 * \tparam t The type to recurse on; must have `std::ranges::value_type_t<rng_t>`.
 *
 * \details
 *
 * Attention, this transformation trait implicitly removes cv-qualifiers on all value_types except the one returned.
 */
template <typename t>
//!\cond
    requires detail::has_value_type<t>
//!\endcond
struct innermost_value_type
{
    //!\brief The return type (recursion not shown).
    using type = value_type_t<remove_cvref_t<t>>;
};

//!\cond
template <typename t>
    requires detail::has_value_type<t> && detail::has_value_type<value_type_t<remove_cvref_t<t>>>
struct innermost_value_type<t>
{
    using type = typename innermost_value_type<value_type_t<remove_cvref_t<t>>>::type;
};
//!\endcond

//!\brief Shortcut for seqan3::innermost_value_type (transformation_trait shortcut).
//!\see seqan3::innermost_value_type
template <typename t>
using innermost_value_type_t = typename innermost_value_type<t>::type;

// ----------------------------------------------------------------------------
// dimension_v
// ----------------------------------------------------------------------------

//NOTE(h-2): this could be moved to a separate file, because it also applies to iterators

/*!\brief Returns the number of times you can call `seqan3::value_type_t` recursively on t (type trait).
 * \tparam t The type to be queried; must resolve `seqan3::value_type_t` at least once.
 *
 * \details
 *
 * Attention, this type trait implicitly removes cv-qualifiers and reference from the types it recurses on and
 * returns.
 */
template <typename t>
//!\cond
    requires detail::has_value_type<t>
//!\endcond
constexpr size_t dimension_v = 1;

//!\cond
template <typename t>
    requires detail::has_value_type<t> && detail::has_value_type<value_type_t<remove_cvref_t<t>>>
constexpr size_t dimension_v<t> = dimension_v<value_type_t<remove_cvref_t<t>>> + 1;
//!\endcond

// ----------------------------------------------------------------------------
// compatible
// ----------------------------------------------------------------------------

//NOTE(h-2): this could be moved to a separate file, because it also applies to iterators

/*!\interface seqan3::compatible <>
 * \brief Two types are "compatible" if their seqan3::dimension_v and their seqan3::innermost_value_type_t are
 * the same.
 *
 * \details
 *
 * \include test/snippet/core/type_traits/range.cpp
 *
 * Attention, this concept implicitly removes cv-qualifiers and reference from the types it recurses on and
 * compares.
 */
//!\cond
template <typename t1, typename t2>
SEQAN3_CONCEPT compatible = requires (t1, t2)
{
    requires (dimension_v<t1> == dimension_v<t2>);

    requires std::is_same_v<innermost_value_type_t<t1>, innermost_value_type_t<t2>>;
};
//!\endcond

//!\}

} // namespace seqan3
