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

#include <seqan3/core/type_list/type_list.hpp>
#include <seqan3/core/type_traits/basic.hpp>

// ----------------------------------------------------------------------------
// seqan3::pack_traits::detail
// ----------------------------------------------------------------------------

namespace seqan3::pack_traits::detail
{

/*!\brief Implementation for seqan3::pack_traits::find.
 * \tparam query_t  The type you are searching for.
 * \tparam pack_t   The type pack.
 * \returns The position of the first occurence of `query_t` in `pack_t` or `-1` if it is not contained.
 * \ingroup type_list
 */
template <typename query_t, typename ...pack_t>
constexpr ptrdiff_t find()
{
    ptrdiff_t i = 0;
    return ((SEQAN3_IS_SAME(query_t, pack_t) ? false : ++i) && ...) ? -1 : i;
}

/*!\brief Implementation for seqan3::pack_traits::find_if.
 * \tparam pred_t   The predicate that is being evaluated.
 * \tparam pack_t   The type pack.
 * \returns The position of the first type `t` in `pack_t` for whom `pred_t<t>::value` is true.
 * \ingroup type_list
 */
template <template <typename> typename pred_t, typename ...pack_t>
constexpr ptrdiff_t find_if()
{
    ptrdiff_t i = 0;
    return ((pred_t<pack_t>::value ? false : ++i) && ...) ? -1 : i;
}

/*!\brief Implementation for seqan3::pack_traits::at.
 * \tparam idx      The index.
 * \tparam head_t   Currently viewed pack_t element.
 * \tparam tail_t   Rest of the type pack.
 * \ingroup type_list
 */
template <ptrdiff_t idx, typename head_t, typename ...tail_t>
auto at()
{
    if constexpr (idx == 0)
        return std::type_identity<head_t>{};
    else if constexpr (idx > 0)
#ifdef __clang__
        return std::type_identity<__type_pack_element<idx - 1, tail_t...>>{};
#else
        return at<idx - 1, tail_t...>();
#endif // __clang__
    else
        return at<sizeof...(tail_t) + 1 + idx, head_t, tail_t...>();
}

/*!\brief Implementation for seqan3::pack_traits::front.
 * \tparam head_t   Currently viewed pack_t element.
 * \tparam tail_t   Rest of the type pack.
 * \ingroup type_list
 */
template <typename head_t, typename ...tail_t>
std::type_identity<head_t> front();

/*!\brief Implementation for seqan3::pack_traits::drop_front.
 * \tparam head_t   Currently viewed pack_t element.
 * \tparam tail_t   Rest of the type pack.
 * \ingroup type_list
 */
template <typename head_t, typename ...tail_t>
type_list<tail_t...> drop_front();

} // namespace seqan3::pack_traits::detail

// ----------------------------------------------------------------------------
// seqan3::pack_traits
// ----------------------------------------------------------------------------

namespace seqan3::pack_traits
{

/*!\name Type pack traits (return a value)
 * \{
 */

/*!\brief The size of a type pack.
 * \tparam pack_t The type pack.
 * \returns `sizeof...(pack_t)`
 * \ingroup type_list
 *
 * \details
 *
 * ### (Compile-time) Complexity
 *
 * * Number of template instantiations: O(1)
 * * Other operations: O(1)
 *
 * \include test/snippet/core/type_list/pack_traits_size.cpp
 */
template <typename ...pack_t>
inline constexpr size_t size = sizeof...(pack_t);

/*!\brief Count the occurrences of a type in a pack.
 * \tparam query_t  The type you are searching for.
 * \tparam pack_t   The type pack.
 * \returns The number of occurrences of the `query_t` in `pack_t`.
 * \ingroup type_list
 *
 * \details
 *
 * ### (Compile-time) Complexity
 *
 * * Number of template instantiations: O(1)
 * * Other operations: O(n)
 *
 * \include test/snippet/core/type_list/pack_traits_count.cpp
 */
template <typename query_t, typename ...pack_t>
inline constexpr ptrdiff_t count =  (SEQAN3_IS_SAME(query_t, pack_t) + ... + 0);

/*!\brief Get the index of the first occurrence of a type in a pack.
 * \tparam query_t  The type you are searching for.
 * \tparam pack_t   The type pack.
 * \returns The position of the first occurrence of `query_t` in `pack_t` or `-1` if it is not contained.
 * \ingroup type_list
 *
 * \details
 *
 * ### (Compile-time) Complexity
 *
 * * Number of template instantiations: O(1)
 * * Other operations: O(n), possibly == `i`, where `i` is the return value
 *
 * \include test/snippet/core/type_list/pack_traits_find.cpp
 */
template <typename query_t, typename ...pack_t>
inline constexpr ptrdiff_t find = seqan3::pack_traits::detail::find<query_t, pack_t...>();

/*!\brief Get the index of the first type in a pack that satisfies the given predicate.
 * \tparam pred_t   The predicate that is being evaluated (a class template).
 * \tparam pack_t   The type pack.
 * \returns The index or `-1` if no types match.
 * \ingroup type_list
 *
 * \details
 *
 * Note that the predicate must be given as a type template (variable templates cannot be passed as template arguments
 * unfortunately). This means e.g. `find_if<std::is_integral, float, double, int, float>` (not `std::is_integral_v`!).
 *
 * ### (Compile-time) Complexity
 *
 * * Number of template instantiations: O(n), possibly == `i`, where `i` is the return value
 * * Other operations: O(n), possibly == `i`, where `i` is the return value
 *
 *  Only the predicate is instantiated.
 *
 * \include test/snippet/core/type_list/pack_traits_find_if.cpp
 */
template <template <typename> typename pred_t, typename ...pack_t>
inline constexpr ptrdiff_t find_if = seqan3::pack_traits::detail::find_if<pred_t, pack_t...>();

/*!\brief Whether a type occurs in a pack or not.
 * \tparam query_t  The type you are searching for.
 * \tparam pack_t   The type pack.
 * \returns `true` or `false`.
 * \ingroup type_list
 *
 * \details
 *
 * ### (Compile-time) Complexity
 *
 * * Number of template instantiations: O(1)
 * * Other operations: O(n), possibly == `i`, where `i` is the index of the first occurrence
 *
 * \include test/snippet/core/type_list/pack_traits_find.cpp
 */
template <typename query_t, typename ...pack_t>
inline constexpr bool contains = (find<query_t, pack_t...> != -1);
//!\}

/*!\name Type pack traits (return a single type)
 * \{
 */

/*!\brief Return the type at given index from the type pack.
 * \tparam idx    The index; must be smaller than the size of the pack.
 * \tparam pack_t The type pack.
 * \ingroup type_list
 *
 * \details
 *
 * Negative indexes are supported (e.g. `at<-1, int, double, bool &>` is `bool &`).
 *
 * ### (Compile-time) Complexity
 *
 * * Number of template instantiations: O(n)
 * * Other operations: O(n)
 *
 * \include test/snippet/core/type_list/pack_traits_at.cpp
 */
template <ptrdiff_t idx, typename ...pack_t>
//!\cond
    requires (idx >= 0 && idx < sizeof...(pack_t)) ||
             (-idx <= sizeof...(pack_t))
//!\endcond
using at = typename decltype(detail::at<idx, pack_t...>())::type;

/*!\brief Return the first type from the type pack.
 * \tparam pack_t The type pack.
 * \ingroup type_list
 *
 * \details
 *
 * ### (Compile-time) Complexity
 *
 * * Number of template instantiations: O(1)
 * * Other operations: O(1)
 *
 * \include test/snippet/core/type_list/pack_traits_front.cpp
 */
template <typename ...pack_t>
//!\cond
    requires (sizeof...(pack_t) > 0)
//!\endcond
using front = typename decltype(detail::front<pack_t...>())::type;

/*!\brief Return the last type from the type pack.
 * \tparam pack_t The type pack.
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
 * \include test/snippet/core/type_list/pack_traits_back.cpp
 */
template <typename ...pack_t>
//!\cond
    requires (sizeof...(pack_t) > 0)
//!\endcond
using back = typename decltype((std::type_identity<pack_t>{}, ...))::type; // use comma operator

//!\}

/*!\name Type pack traits (return a type list)
 * \{
 */

/*!\brief Return a seqan3::type_list of all the types in the type pack, except the first.
 * \tparam pack_t The type pack.
 * \ingroup type_list
 *
 * \details
 *
 * ### (Compile-time) Complexity
 *
 * * Number of template instantiations: O(1)
 * * Other operations: O(1)
 *
 * \include test/snippet/core/type_list/pack_traits_drop_front.cpp
 */
template <typename ...pack_t>
//!\cond
    requires (sizeof...(pack_t) > 0)
//!\endcond
using drop_front = typename decltype(detail::drop_front<pack_t...>())::type;

/*!\brief Apply a transformation trait to every type in the pack and return a seqan3::type_list of the results.
 * \tparam trait_t The transformation trait, **can be an alias template**, e.g. a transformation trait shortcut.
 * \tparam pack_t  The type pack.
 * \ingroup type_list
 *
 * \details
 *
 * The transformation trait given as first argument can be an alias template, e.g. std::type_identity_t, not
 * std::type_identity. The alias must take exactly one argument and be defined for all types in the pack.
 *
 * ### (Compile-time) Complexity
 *
 * * Number of template instantiations: O(n)
 * * Other operations: O(1)
 *
 * \include test/snippet/core/type_list/pack_traits_transform.cpp
 */
template <template <typename> typename trait_t, typename ...pack_t>
using transform = seqan3::type_list<trait_t<pack_t>...>;

//!\}

} // namespace seqan3::pack_traits

// ----------------------------------------------------------------------------
// seqan3::list_traits
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
seqan3::pack_traits::drop_front<pack_t...> drop_front(type_list<pack_t...>);

/*!\brief Implementation for take, drop, take_last and drop_last in seqan3::list_traits.
 * \tparam idx       The index after which to split.
 * \tparam pack1_t   Types in the first type list.
 * \tparam head_t    Head of the second type list.
 * \tparam pack2_t   Rest of the second type list.
 * \ingroup type_list
 */
template <ptrdiff_t idx,
          typename ...pack1_t,
          typename head_t, typename ...pack2_t>
auto split_after(type_list<pack1_t...>, type_list<head_t, pack2_t...>)
{
    if constexpr (idx == sizeof...(pack2_t) + 1)
        return std::pair<type_list<pack1_t..., head_t, pack2_t...>, type_list<>>{};
    else if constexpr (idx == 0)
        return std::pair<type_list<pack1_t...>, type_list<head_t, pack2_t...>>{};
    else
        return split_after<idx - 1>(type_list<pack1_t..., head_t>{}, type_list<pack2_t...>{});
}

/*!\brief Implementation for transform.
 * \tparam trait_t The trait to transform, **must be an alias template**, e.g. a transformation trait shortcut.
 * \tparam pack_t  The type pack.
 * \ingroup type_list
 */
template <template <typename> typename trait_t, typename ...pack_t>
type_list<trait_t<pack_t>...> transform(type_list<pack_t...>);

/*!\brief Implementation for replace_at.
 * \tparam replace_t The type replacing the old one.
 * \tparam pack1_t   Types before the replacement.
 * \tparam head_t    The type to be replaced.
 * \tparam pack2_t   Types after the replacement.
 * \ingroup type_list
 */
template <typename replace_t,
          typename ...pack1_t,
          typename head_t, typename ...pack2_t>
type_list<pack1_t..., replace_t, pack2_t...>
replace_at(std::pair<type_list<pack1_t...>, type_list<head_t, pack2_t...>>);

} // namespace seqan3::list_traits::detail

// ----------------------------------------------------------------------------
// seqan3::list_traits
// ----------------------------------------------------------------------------

//TODO think about whether the type_list_specialisation concept check is expensive and necessary

//!\brief Namespace containing traits for working on seqan3::type_list.
namespace seqan3::list_traits
{

/*!\name Type list traits (return a value)
 * \{
 */

//!\cond
template <seqan3::detail::type_list_specialisation list_t>
inline constexpr size_t size = 0;
//!\endcond

/*!\brief The size of a type list.
 * \ingroup type_list
 * \copydetails seqan3::pack_traits::size
 *
 * \include test/snippet/core/type_list/list_traits_size.cpp
 */
template <typename ...pack_t>
inline constexpr size_t size<type_list<pack_t...>> = sizeof...(pack_t);

//!\cond
template <typename query_t, seqan3::detail::type_list_specialisation list_t>
inline constexpr ptrdiff_t count = -1;
//!\endcond

/*!\brief Count the occurrences of a type in a type list.
 * \ingroup type_list
 * \copydetails seqan3::pack_traits::count
 *
 * \include test/snippet/core/type_list/list_traits_count.cpp
 */
template <typename query_t, typename ...pack_t>
inline constexpr ptrdiff_t count<query_t, type_list<pack_t...>> =
    seqan3::pack_traits::count<query_t, pack_t...>;

//!\cond
template <typename query_t, seqan3::detail::type_list_specialisation list_t>
inline constexpr ptrdiff_t find = -1;
//!\endcond

/*!\brief Get the index of the first occurrence of a type in a type list.
 * \ingroup type_list
 * \copydetails seqan3::pack_traits::find
 *
 * \include test/snippet/core/type_list/list_traits_find.cpp
 */
template <typename query_t, typename ...pack_t>
inline constexpr ptrdiff_t find<query_t, type_list<pack_t...>> =
    seqan3::pack_traits::detail::find<query_t, pack_t...>();

//!\cond
template <template <typename> typename pred_t, seqan3::detail::type_list_specialisation list_t>
inline constexpr ptrdiff_t find_if = -1;
//!\endcond

/*!\brief Get the index of the first type in a type list that satisfies the given predicate.
 * \ingroup type_list
 * \copydetails seqan3::pack_traits::find_if
 *
 * \include test/snippet/core/type_list/list_traits_find.cpp
 */
template <template <typename> typename pred_t, typename ...pack_t>
inline constexpr ptrdiff_t find_if<pred_t, type_list<pack_t...>> =
    seqan3::pack_traits::detail::find_if<pred_t, pack_t...>();

/*!\brief Whether a type occurs in a type list or not.
 * \ingroup type_list
 * \copydetails seqan3::pack_traits::contains
 *
 * \include test/snippet/core/type_list/list_traits_contains.cpp
 */
template <typename query_t, seqan3::detail::type_list_specialisation list_t>
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
 * \include test/snippet/core/type_list/list_traits_at.cpp
 */
template <ptrdiff_t idx, seqan3::detail::type_list_specialisation list_t>
//!\cond
    requires (idx >= 0 && idx < size<list_t>) || (-idx <= size<list_t>)
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
 * \include test/snippet/core/type_list/list_traits_front.cpp
 */
template <seqan3::detail::type_list_specialisation list_t>
//!\cond
    requires (size<list_t> > 0)
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
 * \include test/snippet/core/type_list/list_traits_back.cpp
 */
template <seqan3::detail::type_list_specialisation list_t>
//!\cond
    requires (size<list_t> > 0)
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
 * \include test/snippet/core/type_list/list_traits_concat.cpp
 */
template <typename ...lists_t>
//!\cond
    requires (seqan3::detail::type_list_specialisation<lists_t> && ...)
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
 * \include test/snippet/core/type_list/list_traits_drop_front.cpp
 */
template <seqan3::detail::type_list_specialisation list_t>
//!\cond
    requires (size<list_t> > 0)
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
 * \include test/snippet/core/type_list/list_traits_take.cpp
 */
template <ptrdiff_t i, seqan3::detail::type_list_specialisation list_t>
//!\cond
    requires (i >= 0 && i <= size<list_t>)
//!\endcond
using take = typename decltype(detail::split_after<i>(type_list<>{}, list_t{}))::first_type;

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
 * \include test/snippet/core/type_list/list_traits_drop.cpp
 */
template <ptrdiff_t i, seqan3::detail::type_list_specialisation list_t>
//!\cond
    requires (i >= 0 && i <= size<list_t>)
//!\endcond
using drop = typename decltype(detail::split_after<i>(type_list<>{}, list_t{}))::second_type;

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
 * \include test/snippet/core/type_list/list_traits_take_last.cpp
 */
template <ptrdiff_t i, seqan3::detail::type_list_specialisation list_t>
//!\cond
    requires (i >= 0 && i <= size<list_t>)
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
 * \include test/snippet/core/type_list/list_traits_drop_last.cpp
 */
template <ptrdiff_t i, seqan3::detail::type_list_specialisation list_t>
//!\cond
    requires (i >= 0 && i <= size<list_t>)
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
 * \include test/snippet/core/type_list/list_traits_drop_last.cpp
 */
template <ptrdiff_t i, seqan3::detail::type_list_specialisation list_t>
//!\cond
    requires (i >= 0 && i <= size<list_t>)
//!\endcond
using split_after = decltype(detail::split_after<i>(type_list<>{}, list_t{}));

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
 * \include test/snippet/core/type_list/list_traits_transform.cpp
 */
template <template <typename> typename trait_t, seqan3::detail::type_list_specialisation list_t>
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
 * \include test/snippet/core/type_list/list_traits_replace_at.cpp
 */
template <typename replace_t, std::ptrdiff_t i, seqan3::detail::type_list_specialisation list_t>
//!\cond
    requires (i >= 0 && i < size<list_t>)
//!\endcond
using replace_at = decltype(detail::replace_at<replace_t>(detail::split_after<i>(type_list<>{}, list_t{})));

//!\}

} // namespace seqan3::list_traits

// ----------------------------------------------------------------------------
// seqan3::pack_traits that depend on list_traits
// ----------------------------------------------------------------------------

namespace seqan3::pack_traits
{

/*!\name Type pack traits (return a type list)
 * \{
 */

/*!\brief Return a seqan3::type_list of the first `n` types in the type pack.
 * \tparam i        The target size; must be >= 0 and <= the size of the type pack.
 * \tparam pack_t   The type pack.
 * \ingroup type_list
 *
 * \details
 *
 * ### (Compile-time) Complexity
 *
 * * Number of template instantiations: O(n)
 * * Other operations: O(n)
 *
 * \include test/snippet/core/type_list/pack_traits_take.cpp
 */
template <ptrdiff_t i, typename ...pack_t>
//!\cond
    requires (i >= 0 && i <= sizeof...(pack_t))
//!\endcond
using take = seqan3::list_traits::take<i, type_list<pack_t...>>;

/*!\brief Return a seqan3::type_list of the types in the type pack, except the first `n`.
 * \tparam i        The amount to drop; must be >= 0 and <= the size of the type pack.
 * \tparam pack_t   The type pack.
 * \ingroup type_list
 *
 * \details
 *
 * ### (Compile-time) Complexity
 *
 * * Number of template instantiations: O(n)
 * * Other operations: O(n)
 *
 * \include test/snippet/core/type_list/pack_traits_drop.cpp
 */
template <ptrdiff_t i, typename ...pack_t>
//!\cond
    requires (i >= 0 && i <= sizeof...(pack_t))
//!\endcond
using drop = seqan3::list_traits::drop<i, type_list<pack_t...>>;

/*!\brief Return a seqan3::type_list of the last `n` types in the type pack.
 * \tparam i        The target size; must be >= 0 and <= the size of the type pack.
 * \tparam pack_t   The type pack.
 * \ingroup type_list
 *
 * \details
 *
 * ### (Compile-time) Complexity
 *
 * * Number of template instantiations: O(n)
 * * Other operations: O(n)
 *
 * \include test/snippet/core/type_list/pack_traits_take_last.cpp
 */
template <ptrdiff_t i, typename ...pack_t>
//!\cond
    requires (i >= 0 && i <= sizeof...(pack_t))
//!\endcond
using take_last = seqan3::list_traits::take_last<i, type_list<pack_t...>>;

/*!\brief Return a seqan3::type_list of the types the type pack, except the last `n`.
 * \tparam i        The amount to drop; must be >= 0 and <= the size of the type pack.
 * \tparam pack_t   The type pack.
 * \ingroup type_list
 *
 * \details
 *
 * ### (Compile-time) Complexity
 *
 * * Number of template instantiations: O(n)
 * * Other operations: O(n)
 *
 * \include test/snippet/core/type_list/pack_traits_take_last.cpp
 */
template <ptrdiff_t i, typename ...pack_t>
//!\cond
    requires (i >= 0 && i <= sizeof...(pack_t))
//!\endcond
using drop_last = seqan3::list_traits::drop_last<i, type_list<pack_t...>>;

/*!\brief Split a type pack into two parts returned as a pair of seqan3::type_list.
 * \tparam i        The number of elements after which to split; must be >= 0 and <= the size of the type pack.
 * \tparam pack_t   The type pack.
 * \ingroup type_list
 *
 * \details
 *
 * ### (Compile-time) Complexity
 *
 * * Number of template instantiations: O(n)
 * * Other operations: O(n)
 *
 * \include test/snippet/core/type_list/pack_traits_take_last.cpp
 */
template <ptrdiff_t i, typename ...pack_t>
//!\cond
    requires (i >= 0 && i <= sizeof...(pack_t))
//!\endcond
using split_after = seqan3::list_traits::split_after<i, type_list<pack_t...>>;

/*!\brief Replace the type at the given index with the given type.
 * \tparam replace_t The type to replace the old type with.
 * \tparam i         The index of the type to be replaced.
 * \tparam pack_t    The (input) type pack.
 * \ingroup type_list
 *
 * \details
 *
 * ### (Compile-time) Complexity
 *
 * * Number of template instantiations: O(n)
 * * Other operations: O(n)
 *
 * \include test/snippet/core/type_list/pack_traits_take_last.cpp
 */
template <typename replace_t, std::ptrdiff_t i, typename ...pack_t>
//!\cond
    requires (i >= 0 && i < sizeof...(pack_t))
//!\endcond
using replace_at = decltype(seqan3::list_traits::replace_at<replace_t, i, type_list<pack_t...>>());

//!\}

} // namespace seqan3::pack_traits
