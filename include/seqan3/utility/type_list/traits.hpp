// -----------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2020, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2020, Knut Reinert & MPI für molekulare Genetik
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at: https://github.com/seqan/seqan3/blob/master/LICENSE.md
// -----------------------------------------------------------------------------------------------------

/*!\file
 * \author Hannes Hauswedell <hannes.hauswedell AT fu-berlin.de>
 * \brief Provides traits for seqan3::type_list.
 */

#pragma once

#include <seqan3/std/type_traits>

#include <seqan3/utility/type_list/type_list.hpp>
#include <seqan3/utility/type_pack/traits.hpp>

// ----------------------------------------------------------------------------
// seqan3::list_traits::detail
// ----------------------------------------------------------------------------

namespace seqan3::list_traits::detail
{

/*!\brief Implementation for seqan3::list_traits::at.
 * \tparam idx      The index.
 * \tparam pack_t   Types in the type list.
 * \ingroup type_list
 */
template <ptrdiff_t idx, typename ...pack_t>
std::type_identity<seqan3::pack_traits::at<idx, pack_t...>> at(type_list<pack_t...>);

/*!\brief Implementation for seqan3::list_traits::front.
 * \tparam pack_t   Types in the type list.
 * \ingroup type_list
 */
template <typename ...pack_t>
std::type_identity<seqan3::pack_traits::front<pack_t...>> front(type_list<pack_t...>);

/*!\brief Implementation for seqan3::list_traits::back.
 * \tparam pack_t   Types in the type list.
 * \ingroup type_list
 */
template <typename ...pack_t>
std::type_identity<seqan3::pack_traits::back<pack_t...>> back(type_list<pack_t...>);

/*!\brief Implementation for seqan3::list_traits::concat.
 * \tparam pack1_t   Types in the first type list.
 * \tparam pack2_t   Types in the second type list.
 * \ingroup type_list
 */
template <typename ...pack1_t,
          typename ...pack2_t>
type_list<pack1_t..., pack2_t...> concat(type_list<pack1_t...>, type_list<pack2_t...>);

/*!\brief Implementation for seqan3::list_traits::concat [overload for more than two lists].
 * \tparam pack1_t      Types in the first type list.
 * \tparam pack2_t      Types in the second type list.
 * \tparam more_lists_t The remaining type lists.
 * \ingroup type_list
 */
template <typename ...pack1_t,
          typename ...pack2_t,
          typename ...more_lists_t>
auto concat(type_list<pack1_t...>, type_list<pack2_t...>, more_lists_t ...)
{
    return concat(type_list<pack1_t..., pack2_t...>{}, more_lists_t{}...);
}

/*!\brief Implementation for seqan3::list_traits::drop_front.
 * \tparam pack_t   Types in the type list.
 * \ingroup type_list
 */
template <typename ...pack_t>
pack_traits::drop_front<pack_t...> drop_front(type_list<pack_t...>);

/*!\brief Implementation for seqan3::list_traits::transform.
 * \tparam trait_t The trait to transform, **must be an alias template**, e.g. a transformation trait shortcut.
 * \tparam pack_t  Types in the type list.
 * \ingroup type_list
 */
template <template <typename> typename trait_t, typename ...pack_t>
pack_traits::transform<trait_t, pack_t...> transform(type_list<pack_t...>);

/*!\brief Implementation for seqan3::list_traits::split_after.
 * \tparam idx The index after which to split.
 * \tparam pack_t Types in the type list to split
 * \ingroup type_list
 */
template <ptrdiff_t idx,
          typename ...pack1_t>
pack_traits::split_after<idx, pack1_t...> split_after(type_list<pack1_t...>);

/*!\brief Implementation for seqan3::list_traits::replace_at.
 * \tparam replace_t The type replacing the old one.
 * \tparam idx The index of the type to replace.
 * \tparam pack_t Types in the type list to be modified.
 * \ingroup type_list
 */
template <typename replace_t,
          ptrdiff_t idx,
          typename ...pack_t>
pack_traits::replace_at<replace_t, idx, pack_t...> replace_at(type_list<pack_t...>);

//!\brief A replacement for meta::reverse [recursion anchor]
inline constexpr type_list<> reverse(type_list<>) { return {}; }

//!\brief A replacement for meta::reverse [recursion]
template <typename head_t, typename ...pack_t>
auto reverse(type_list<head_t, pack_t...>)
{
    return concat(reverse(type_list<pack_t...>{}), type_list<head_t>{});
}

} // namespace seqan3::list_traits::detail

// ----------------------------------------------------------------------------
// seqan3::list_traits
// ----------------------------------------------------------------------------

//!\brief Namespace containing traits for working on seqan3::type_list.
namespace seqan3::list_traits
{

/*!\name Type list traits (return a value)
 * \{
 */

//!\cond
template <typename list_t>
    requires seqan3::detail::template_specialisation_of<list_t, seqan3::type_list>
inline constexpr size_t size = 0;
//!\endcond

/*!\brief The size of a type list.
 * \ingroup type_list
 * \copydetails seqan3::pack_traits::size
 *
 * \include test/snippet/utility/type_list/list_traits_size.cpp
 */
template <typename ...pack_t>
inline constexpr size_t size<type_list<pack_t...>> = sizeof...(pack_t);

//!\cond
template <typename query_t, typename list_t>
    requires seqan3::detail::template_specialisation_of<list_t, seqan3::type_list>
inline constexpr ptrdiff_t count = -1;
//!\endcond

/*!\brief Count the occurrences of a type in a type list.
 * \ingroup type_list
 * \copydetails seqan3::pack_traits::count
 *
 * \include test/snippet/utility/type_list/list_traits_count.cpp
 */
template <typename query_t, typename ...pack_t>
inline constexpr ptrdiff_t count<query_t, type_list<pack_t...>> =
    seqan3::pack_traits::count<query_t, pack_t...>;

//!\cond
template <typename query_t, typename list_t>
    requires seqan3::detail::template_specialisation_of<list_t, seqan3::type_list>
inline constexpr ptrdiff_t find = -1;
//!\endcond

/*!\brief Get the index of the first occurrence of a type in a type list.
 * \ingroup type_list
 * \copydetails seqan3::pack_traits::find
 *
 * \include test/snippet/utility/type_list/list_traits_find.cpp
 */
template <typename query_t, typename ...pack_t>
inline constexpr ptrdiff_t find<query_t, type_list<pack_t...>> =
    seqan3::pack_traits::detail::find<query_t, pack_t...>();

//!\cond
template <template <typename> typename pred_t, typename list_t>
    requires seqan3::detail::template_specialisation_of<list_t, seqan3::type_list>
inline constexpr ptrdiff_t find_if = -1;
//!\endcond

/*!\brief Get the index of the first type in a type list that satisfies the given predicate.
 * \ingroup type_list
 * \copydetails seqan3::pack_traits::find_if
 *
 * \include test/snippet/utility/type_list/list_traits_find.cpp
 */
template <template <typename> typename pred_t, typename ...pack_t>
inline constexpr ptrdiff_t find_if<pred_t, type_list<pack_t...>> =
    seqan3::pack_traits::detail::find_if<pred_t, pack_t...>();

/*!\brief Whether a type occurs in a type list or not.
 * \ingroup type_list
 * \copydetails seqan3::pack_traits::contains
 *
 * \include test/snippet/utility/type_list/list_traits_contains.cpp
 */
template <typename query_t, typename list_t>
//!\cond
    requires seqan3::detail::template_specialisation_of<list_t, seqan3::type_list>
//!\endcond
inline constexpr bool contains = (find<query_t, list_t> != -1);

//!\}

/*!\name Type list traits (return a single type)
 * \{
 */

/*!\brief Return the type at given index from the type list.
 * \tparam idx    The index; must be smaller than the size of the type list.
 * \tparam list_t The type_list.
 * \ingroup type_list
 *
 * \details
 *
 * Negative indexes are supported (e.g. `at<-1, type_list<int, double, bool &>>` is `bool &`).
 *
 * ### (Compile-time) Complexity
 *
 * * Number of template instantiations: O(n)
 * * Other operations: O(n)
 *
 * \include test/snippet/utility/type_list/list_traits_at.cpp
 */
template <ptrdiff_t idx, typename list_t>
//!\cond
    requires (seqan3::detail::template_specialisation_of<list_t, seqan3::type_list>) &&
             ((idx >= 0 && idx < size<list_t>) || (-idx <= size<list_t>))
//!\endcond
using at = typename decltype(detail::at<idx>(list_t{}))::type;

/*!\brief Return the first type from the type list.
 * \tparam list_t The type list.
 * \ingroup type_list
 *
 * \details
 *
 * ### (Compile-time) Complexity
 *
 * * Number of template instantiations: O(1)
 * * Other operations: O(1)
 *
 * \include test/snippet/utility/type_list/list_traits_front.cpp
 */
template <typename list_t>
//!\cond
    requires (seqan3::detail::template_specialisation_of<list_t, seqan3::type_list>) && (size<list_t> > 0)
//!\endcond
using front = typename decltype(detail::front(list_t{}))::type;

/*!\brief Return the last type from the type list.
 * \tparam list_t The type list.
 * \ingroup type_list
 *
 * \details
 *
 * ### (Compile-time) Complexity
 *
 * * Number of template instantiations: O(n) (possibly O(1))
 * * Other operations: O(1)
 *
 * Notably faster than `seqan3::pack_traits::at<size<pack...> - 1, pack...>` (no recursive template
 * instantiations).
 *
 * \include test/snippet/utility/type_list/list_traits_back.cpp
 */
template <typename list_t>
//!\cond
    requires (seqan3::detail::template_specialisation_of<list_t, seqan3::type_list>) && (size<list_t> > 0)
//!\endcond
using back = typename decltype(detail::back(list_t{}))::type;

//!\}

/*!\name Type list traits (return a type list)
 * \{
 */

/*!\brief Join two seqan3::type_list s into one.
 * \tparam list1_t The first (input) type list.
 * \tparam list2_t The second (input) type list.
 * \ingroup type_list
 *
 * \details
 *
 * ### (Compile-time) Complexity
 *
 * * Number of template instantiations: O(n) in the number of type lists
 * * Other operations: O(n) in the number of type lists
 *
 * Complexity is independent of the number of types in each list.
 *
 * \include test/snippet/utility/type_list/list_traits_concat.cpp
 */
template <typename ...lists_t>
//!\cond
    requires (seqan3::detail::template_specialisation_of<lists_t, seqan3::type_list> && ...)
//!\endcond
using concat = decltype(detail::concat(lists_t{}...));

/*!\brief Return a seqan3::type_list of all the types in the type list, except the first.
 * \tparam list_t The (input) type list.
 * \ingroup type_list
 *
 * \details
 *
 * ### (Compile-time) Complexity
 *
 * * Number of template instantiations: O(1)
 * * Other operations: O(1)
 *
 * \include test/snippet/utility/type_list/list_traits_drop_front.cpp
 */
template <typename list_t>
//!\cond
    requires (seqan3::detail::template_specialisation_of<list_t, seqan3::type_list>) && (size<list_t> > 0)
//!\endcond
using drop_front = decltype(detail::drop_front(list_t{}));

/*!\brief Return a seqan3::type_list of the first `n` types in the input type list.
 * \tparam i        The target size; must be >= 0 and <= the size of the input type list.
 * \tparam list_t   The (input) type list.
 * \ingroup type_list
 *
 * \details
 *
 * ### (Compile-time) Complexity
 *
 * * Number of template instantiations: O(n)
 * * Other operations: O(n)
 *
 * \include test/snippet/utility/type_list/list_traits_take.cpp
 */
template <ptrdiff_t i, typename list_t>
//!\cond
    requires (seqan3::detail::template_specialisation_of<list_t, seqan3::type_list>) && (i >= 0 && i <= size<list_t>)
//!\endcond
using take = typename decltype(detail::split_after<i>(list_t{}))::first_type;

/*!\brief Return a seqan3::type_list of the types in the input type list, except the first `n`.
 * \tparam i        The amount to drop; must be >= 0 and <= the size of the input type list.
 * \tparam list_t   The (input) type list.
 * \ingroup type_list
 *
 * \details
 *
 * ### (Compile-time) Complexity
 *
 * * Number of template instantiations: O(n)
 * * Other operations: O(n)
 *
 * \include test/snippet/utility/type_list/list_traits_drop.cpp
 */
template <ptrdiff_t i, typename list_t>
//!\cond
    requires (seqan3::detail::template_specialisation_of<list_t, seqan3::type_list>) && (i >= 0 && i <= size<list_t>)
//!\endcond
using drop = typename decltype(detail::split_after<i>(list_t{}))::second_type;

/*!\brief Return a seqan3::type_list of the last `n` types in the input type list.
 * \tparam i        The target size; must be >= 0 and <= the size of the input type list.
 * \tparam list_t   The (input) type list.
 * \ingroup type_list
 *
 * \details
 *
 * ### (Compile-time) Complexity
 *
 * * Number of template instantiations: O(n)
 * * Other operations: O(n)
 *
 * \include test/snippet/utility/type_list/list_traits_take_last.cpp
 */
template <ptrdiff_t i, typename list_t>
//!\cond
    requires (seqan3::detail::template_specialisation_of<list_t, seqan3::type_list>) && (i >= 0 && i <= size<list_t>)
//!\endcond
using take_last = drop<size<list_t> - i, list_t>;

/*!\brief Return a seqan3::type_list of the types the input type list, except the last `n`.
 * \tparam i        The amount to drop; must be >= 0 and <= the size of the input type list.
 * \tparam list_t   The (input) type list.
 * \ingroup type_list
 *
 * \details
 *
 * ### (Compile-time) Complexity
 *
 * * Number of template instantiations: O(n)
 * * Other operations: O(n)
 *
 * \include test/snippet/utility/type_list/list_traits_drop_last.cpp
 */
template <ptrdiff_t i, typename list_t>
//!\cond
    requires (seqan3::detail::template_specialisation_of<list_t, seqan3::type_list>) && (i >= 0 && i <= size<list_t>)
//!\endcond
using drop_last = take<size<list_t> - i, list_t>;

/*!\brief Split a seqan3::type_list into two parts returned as a pair of seqan3::type_list.
 * \tparam i        The number of elements after which to split; must be >= 0 and <= the size of the input type list.
 * \tparam list_t   The (input) type list.
 * \ingroup type_list
 *
 * \details
 *
 * ### (Compile-time) Complexity
 *
 * * Number of template instantiations: O(n)
 * * Other operations: O(n)
 *
 * \include test/snippet/utility/type_list/list_traits_drop_last.cpp
 */
template <ptrdiff_t i, typename list_t>
//!\cond
    requires (seqan3::detail::template_specialisation_of<list_t, seqan3::type_list>) && (i >= 0 && i <= size<list_t>)
//!\endcond
using split_after = decltype(detail::split_after<i>(list_t{}));

/*!\brief Apply a transformation trait to every type in the list and return a seqan3::type_list of the results.
 * \tparam trait_t The trait to transform, **must be an alias template**, e.g. a transformation trait shortcut.
 * \tparam list_t  The (input) type list.
 * \ingroup type_list
 *
 * \details
 *
 * The transformation trait given as first argument must be an alias template, e.g. std::type_identity_t, not
 * std::type_identity. The alias must take exactly one argument and be defined for all types in the input list.
 *
 * ### (Compile-time) Complexity
 *
 * * Number of template instantiations: O(n)
 * * Other operations: O(n)
 *
 * \include test/snippet/utility/type_list/list_traits_transform.cpp
 */
template <template <typename> typename trait_t, typename list_t>
//!\cond
    requires seqan3::detail::template_specialisation_of<list_t, seqan3::type_list>
//!\endcond
using transform = decltype(detail::transform<trait_t>(list_t{}));

/*!\brief Replace the type at the given index with the given type.
 * \tparam replace_t The type to replace the old type with.
 * \tparam i         The index of the type to be replaced.
 * \tparam list_t    The (input) type list.
 * \ingroup type_list
 *
 * \details
 *
 * ### (Compile-time) Complexity
 *
 * * Number of template instantiations: O(n)
 * * Other operations: O(n)
 *
 * \include test/snippet/utility/type_list/list_traits_replace_at.cpp
 */
template <typename replace_t, std::ptrdiff_t i, typename list_t>
//!\cond
    requires (seqan3::detail::template_specialisation_of<list_t, seqan3::type_list>) && (i >= 0 && i < size<list_t>)
//!\endcond
using replace_at = decltype(detail::replace_at<replace_t, i>(list_t{}));

//!\}

} // namespace seqan3::list_traits
