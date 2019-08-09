// -----------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2019, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2019, Knut Reinert & MPI für molekulare Genetik
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at: https://github.com/seqan/seqan3/blob/master/LICENSE.md
// -----------------------------------------------------------------------------------------------------

/*!\file
 * \author Hannes Hauswedell <hannes.hauswedell AT fu-berlin.de>
 * \brief Provides traits for seqan3::type_list.
 */

#pragma once

#include <seqan3/core/type_list/type_list.hpp>

// ----------------------------------------------------------------------------
// seqan3::pack_traits::detail
// ----------------------------------------------------------------------------

namespace seqan3::pack_traits::detail
{

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
#if __clang__
        return std::type_identity<__type_pack_element<idx - 1, tail_t...>>;
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
 * Constant
 */
template <typename ...pack_t>
inline constexpr size_t size = sizeof...(pack_t);

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
 * Equal to `idx` (linear).
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
 * Constant.
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
 * Constant.
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
 * Constant.
 */
template <typename ...pack_t>
//!\cond
    requires (sizeof...(pack_t) > 0)
//!\endcond
using drop_front = typename decltype(detail::drop_front<pack_t...>())::type;

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

} // namespace seqan3::list_traits::detail

// ----------------------------------------------------------------------------
// seqan3::list_traits
// ----------------------------------------------------------------------------

//TODO think about whether the TypeList concept check is expensive and necessary

//!\brief Namespace containing traits for working on seqan3::type_list.
namespace seqan3::list_traits
{

//!\cond
template <seqan3::detail::TypeList list_t>
inline constexpr size_t size = 0;
//!\endcond

/*!\name Type list traits (return a value)
 * \{
 */

/*!\brief The size of a type list.
 * \tparam pack_t Types in the list.
 * \returns `sizeof...(pack_t)`
 * \ingroup type_list
 *
 * \details
 *
 * ### (Compile-time) Complexity
 *
 * Constant
 */
template <typename ...pack_t>
inline constexpr size_t size<type_list<pack_t...>> = sizeof...(pack_t);

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
 * Equal to `idx` (linear).
 */
template <ptrdiff_t idx, seqan3::detail::TypeList list_t>
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
 * Constant.
 */
template <seqan3::detail::TypeList list_t>
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
 * Constant.
 */
template <seqan3::detail::TypeList list_t>
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
 * Constant.
 */
template <seqan3::detail::TypeList list1_t, seqan3::detail::TypeList list2_t>
using concat = decltype(detail::concat(list1_t{}, list2_t{}));

/*!\brief Return a seqan3::type_list of all the types in the type list, except the first.
 * \tparam list_t The (input) type list.
 * \ingroup type_list
 *
 * \details
 *
 * ### (Compile-time) Complexity
 *
 * Constant.
 */
template <seqan3::detail::TypeList list_t>
//!\cond
    requires (size<list_t> > 0)
//!\endcond
using drop_front = decltype(detail::drop_front(list_t{}));

/*!\brief Return a seqan3::type_list of the first `n` types in the input type list.
 * \tparam n        The target size; must be >= 0 and <= the size of the input type list.
 * \tparam list_t   The (input) type list.
 * \ingroup type_list
 *
 * \details
 *
 * ### (Compile-time) Complexity
 *
 * Equal to `n` (linear).
 */
template <ptrdiff_t n, seqan3::detail::TypeList list_t>
//!\cond
    requires (n >= 0 && n <= size<list_t>)
//!\endcond
using take = typename decltype(detail::split_after<n>(type_list<>{}, list_t{}))::first_type;

/*!\brief Return a seqan3::type_list of the types in the input type list, except the first `n`.
 * \tparam n        The amount to drop; must be >= 0 and <= the size of the input type list.
 * \tparam list_t   The (input) type list.
 * \ingroup type_list
 *
 * \details
 *
 * ### (Compile-time) Complexity
 *
 * Equal to `n` (linear).
 */
template <ptrdiff_t n, seqan3::detail::TypeList list_t>
//!\cond
    requires (n >= 0 && n <= size<list_t>)
//!\endcond
using drop = typename decltype(detail::split_after<n>(type_list<>{}, list_t{}))::second_type;

/*!\brief Return a seqan3::type_list of the last `n` types in the input type list.
 * \tparam n        The target size; must be >= 0 and <= the size of the input type list.
 * \tparam list_t   The (input) type list.
 * \ingroup type_list
 *
 * \details
 *
 * ### (Compile-time) Complexity
 *
 * Equal to `size<list_t> - n` (linear).
 */
template <ptrdiff_t n, seqan3::detail::TypeList list_t>
//!\cond
    requires (n >= 0 && n <= size<list_t>)
//!\endcond
using take_last = drop<size<list_t> - n, list_t>;

/*!\brief Return a seqan3::type_list of the types the input type list, except the last `n`.
 * \tparam n        The amount to drop; must be >= 0 and <= the size of the input type list.
 * \tparam list_t   The (input) type list.
 * \ingroup type_list
 *
 * \details
 *
 * ### (Compile-time) Complexity
 *
 * Equal to `size<list_t> - n` (linear).
 */
template <ptrdiff_t n, seqan3::detail::TypeList list_t>
//!\cond
    requires (n >= 0 && n <= size<list_t>)
//!\endcond
using drop_last = take<size<list_t> - n, list_t>;

/*!\brief Split a seqan3::type_list into two parts returned as a pair of seqan3::type_list.
 * \tparam n        The number of elements after which to split; must be >= 0 and <= the size of the input type list.
 * \tparam list_t   The (input) type list.
 * \ingroup type_list
 *
 * \details
 *
 * ### (Compile-time) Complexity
 *
 * Equal to `n` (linear).
 */
template <ptrdiff_t n, seqan3::detail::TypeList list_t>
//!\cond
    requires (n >= 0 && n <= size<list_t>)
//!\endcond
using split_after = decltype(detail::split_after<n>(type_list<>{}, list_t{}));

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
 * \tparam n        The target size; must be >= 0 and <= the size of the type pack.
 * \tparam pack_t   The type pack.
 * \ingroup type_list
 *
 * \details
 *
 * ### (Compile-time) Complexity
 *
 * Equal to `n` (linear).
 */
template <ptrdiff_t n, typename ...pack_t>
//!\cond
    requires (n >= 0 && n <= sizeof...(pack_t))
//!\endcond
using take = seqan3::list_traits::take<n, type_list<pack_t...>>;

/*!\brief Return a seqan3::type_list of the types in the type pack, except the first `n`.
 * \tparam n        The amount to drop; must be >= 0 and <= the size of the type pack.
 * \tparam pack_t   The type pack.
 * \ingroup type_list
 *
 * \details
 *
 * ### (Compile-time) Complexity
 *
 * Equal to `n` (linear).
 */
template <ptrdiff_t n, typename ...pack_t>
//!\cond
    requires (n >= 0 && n <= sizeof...(pack_t))
//!\endcond
using drop = seqan3::list_traits::drop<n, type_list<pack_t...>>;

/*!\brief Return a seqan3::type_list of the last `n` types in the type pack.
 * \tparam n        The target size; must be >= 0 and <= the size of the type pack.
 * \tparam pack_t   The type pack.
 * \ingroup type_list
 *
 * \details
 *
 * ### (Compile-time) Complexity
 *
 * Equal to `sizeof...(pack_t) - n` (linear).
 */
template <ptrdiff_t n, typename ...pack_t>
//!\cond
    requires (n >= 0 && n <= sizeof...(pack_t))
//!\endcond
using take_last = seqan3::list_traits::take_last<n, type_list<pack_t...>>;

/*!\brief Return a seqan3::type_list of the types the type pack, except the last `n`.
 * \tparam n        The amount to drop; must be >= 0 and <= the size of the type pack.
 * \tparam pack_t   The type pack.
 * \ingroup type_list
 *
 * \details
 *
 * ### (Compile-time) Complexity
 *
 * Equal to `sizeof...(pack_t) - n (linear).
 */
template <ptrdiff_t n, typename ...pack_t>
//!\cond
    requires (n >= 0 && n <= sizeof...(pack_t))
//!\endcond
using drop_last = seqan3::list_traits::drop_last<n, type_list<pack_t...>>;

/*!\brief Split a type pack into two parts returned as a pair of seqan3::type_list.
 * \tparam n        The number of elements after which to split; must be >= 0 and <= the size of the type pack.
 * \tparam pack_t   The type pack.
 * \ingroup type_list
 *
 * \details
 *
 * ### (Compile-time) Complexity
 *
 * Equal to `n` (linear).
 */
template <ptrdiff_t n, typename ...pack_t>
//!\cond
    requires (n >= 0 && n <= sizeof...(pack_t))
//!\endcond
using split_after = seqan3::list_traits::split_after<n, type_list<pack_t...>>;

//!\}

} // namespace seqan3::pack_traits
