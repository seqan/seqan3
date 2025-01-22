// SPDX-FileCopyrightText: 2006-2025 Knut Reinert & Freie Universität Berlin
// SPDX-FileCopyrightText: 2016-2025 Knut Reinert & MPI für molekulare Genetik
// SPDX-License-Identifier: BSD-3-Clause

/*!\file
 * \brief Provides various transformation traits used by the range module.
 * \author Hannes Hauswedell <hannes.hauswedell AT fu-berlin.de>
 */

#pragma once

#include <iterator>
#include <ranges>
#include <type_traits>

#include <seqan3/core/detail/iterator_traits.hpp>
#include <seqan3/core/platform.hpp>
#include <seqan3/utility/type_traits/basic.hpp>

// TODO(h-2): add innermost_reference instead of or addition to range_innermost_value?

//NOTE(h-2): for the range overloads we explicitly forbid that the type is iteratoer
// because some types are actually both (e.g. std::directory_iterator)

namespace seqan3::detail
{

//!\cond
template <typename t>
concept has_range_value_type = requires { typename std::ranges::range_value_t<std::remove_cvref_t<t>>; };
//!\endcond

//!\brief Makes range_t const if const_range is true; otherwise keeps range_t as is.
//!\ingroup core_range
template <bool const_range, typename range_t>
using maybe_const_range_t = std::conditional_t<const_range, range_t const, range_t>;

//!\brief Returns the const iterator of range_t if const_range is true; otherwise the non-const iterator.
//!\ingroup core_range
template <bool const_range, typename range_t>
using maybe_const_iterator_t = std::ranges::iterator_t<maybe_const_range_t<const_range, range_t>>;

//!\brief Returns the const sentinel of range_t if const_range is true; otherwise the non-const sentinel.
//!\ingroup core_range
template <bool const_v, typename range_t>
using maybe_const_sentinel_t = std::ranges::sentinel_t<maybe_const_range_t<const_v, range_t>>;

//!\cond
//!\brief The same as std::indirect_unary_predicate but on a range type instead of an iterator type.
template <typename unary_predicate_fn_t, typename urng_t>
concept indirect_unary_predicate_on_range =
    std::ranges::range<urng_t> && std::indirect_unary_predicate<unary_predicate_fn_t, std::ranges::iterator_t<urng_t>>;
//!\endcond
} // namespace seqan3::detail

namespace seqan3
{

// ----------------------------------------------------------------------------
// range_innermost_value
// ----------------------------------------------------------------------------

/*!\brief Recursively determines the `value_type` on containers and/or iterators.
 * \ingroup core_range
 * \implements seqan3::transformation_trait
 * \tparam t The type to recurse on; must have `std::ranges::value_type_t<rng_t>`.
 *
 * \details
 *
 * Attention, this transformation trait implicitly removes cv-qualifiers on all value_types except the one returned.
 * \see core_range
 */
template <typename t>
    requires detail::has_range_value_type<t>
struct range_innermost_value
{
    //!\brief The return type (recursion not shown).
    using type = std::ranges::range_value_t<std::remove_cvref_t<t>>;
};

//!\cond
template <typename t>
    requires detail::has_range_value_type<t>
          && detail::has_range_value_type<std::ranges::range_value_t<std::remove_cvref_t<t>>>
struct range_innermost_value<t>
{
    using type = typename range_innermost_value<std::ranges::range_value_t<std::remove_cvref_t<t>>>::type;
};
//!\endcond

//!\brief Shortcut for seqan3::range_innermost_value (transformation_trait shortcut).
//!\ingroup core_range
//!\see seqan3::range_innermost_value
template <typename t>
using range_innermost_value_t = typename range_innermost_value<t>::type;

// ----------------------------------------------------------------------------
// range_dimension_v
// ----------------------------------------------------------------------------

/*!\brief Returns the number of times you can call `seqan3::value_type_t` recursively on t (type trait).
 * \ingroup core_range
 * \tparam t The type to be queried; must resolve `seqan3::value_type_t` at least once.
 *
 * \details
 *
 * Attention, this type trait implicitly removes cv-qualifiers and reference from the types it recurses on and
 * returns.
 */
template <typename t>
    requires detail::has_range_value_type<t>
constexpr size_t range_dimension_v = 1;

//!\cond
template <typename t>
    requires detail::has_range_value_type<t>
              && detail::has_range_value_type<std::ranges::range_value_t<std::remove_cvref_t<t>>>
constexpr size_t range_dimension_v<t> = range_dimension_v<std::ranges::range_value_t<std::remove_cvref_t<t>>> + 1;
//!\endcond

} // namespace seqan3
