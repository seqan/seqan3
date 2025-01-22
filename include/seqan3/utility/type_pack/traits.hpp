// SPDX-FileCopyrightText: 2006-2025 Knut Reinert & Freie Universität Berlin
// SPDX-FileCopyrightText: 2016-2025 Knut Reinert & MPI für molekulare Genetik
// SPDX-License-Identifier: BSD-3-Clause

/*!\file
 * \author Rene Rahn <rene.rahn AT fu-berlin.de>
 * \brief Provides various traits for template packs.
 */

#pragma once

#include <type_traits>
#include <utility>

#include <seqan3/utility/type_list/type_list.hpp>

// ----------------------------------------------------------------------------
// seqan3::pack_traits::detail
// ----------------------------------------------------------------------------

namespace seqan3::pack_traits::detail
{

/*!\brief Implementation for seqan3::pack_traits::find.
 * \tparam query_t  The type you are searching for.
 * \tparam pack_t   The type pack.
 * \returns The position of the first occurrence of `query_t` in `pack_t` or `-1` if it is not contained.
 * \ingroup utility_type_pack
 */
template <typename query_t, typename... pack_t>
constexpr ptrdiff_t find()
{
    ptrdiff_t i = 0;
    return ((std::is_same_v<query_t, pack_t> ? false : ++i) && ...) ? -1 : i;
}

/*!\brief Implementation for seqan3::pack_traits::find_if.
 * \tparam pred_t   The predicate that is being evaluated.
 * \tparam pack_t   The type pack.
 * \returns The position of the first type `t` in `pack_t` for whom ``pred_t<t>::%value`` is true.
 * \ingroup utility_type_pack
 */
template <template <typename> typename pred_t, typename... pack_t>
constexpr ptrdiff_t find_if()
{
    ptrdiff_t i = 0;
    return ((pred_t<pack_t>::value ? false : ++i) && ...) ? -1 : i;
}

/*!\brief Implementation for seqan3::pack_traits::at.
 * \tparam idx      The index.
 * \tparam head_t   Currently viewed pack_t element.
 * \tparam tail_t   Rest of the type pack.
 * \ingroup utility_type_pack
 */
template <ptrdiff_t idx, typename head_t, typename... tail_t>
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
 * \ingroup utility_type_pack
 */
template <typename head_t, typename... tail_t>
std::type_identity<head_t> front();

/*!\brief Implementation for seqan3::pack_traits::drop_front.
 * \tparam head_t   Currently viewed pack_t element.
 * \tparam tail_t   Rest of the type pack.
 * \ingroup utility_type_pack
 */
template <typename head_t, typename... tail_t>
type_list<tail_t...> drop_front();

/*!\brief Implementation for seqan3::pack_traits::split_after.
 * \tparam idx     The index to split the type pack at.
 * \tparam pack1_t The type pack before the split index.
 * \tparam head_t  The next type that is moved before split index.
 * \tparam pack2_t The type pack after the split index.
 * \ingroup utility_type_pack
 */
template <ptrdiff_t idx, typename head_t, typename... pack2_t, typename... pack1_t>
auto split_after(type_list<pack1_t...>)
{
    if constexpr (idx == sizeof...(pack2_t) + 1)
        return std::pair<type_list<pack1_t..., head_t, pack2_t...>, type_list<>>{};
    else if constexpr (idx == 0)
        return std::pair<type_list<pack1_t...>, type_list<head_t, pack2_t...>>{};
    else
        return split_after<idx - 1, pack2_t...>(type_list<pack1_t..., head_t>{});
}

/*!\brief Implementation for seqan3::pack_traits::replace_at.
 * \tparam replace_t The type replacing the old one.
 * \tparam idx The index of the type to replace.
 * \tparam pack_t The type pack to be modified.
 * \tparam i The indicies of the index sequence associated with the type pack.
 * \ingroup utility_type_pack
 */
template <typename replace_t, ptrdiff_t idx, typename... pack_t, size_t... i>
auto replace_at(std::index_sequence<i...>) -> type_list<std::conditional_t<i == idx, replace_t, pack_t>...>;

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
 * \ingroup utility_type_pack
 *
 * \details
 *
 * ### (Compile-time) Complexity
 *
 * * Number of template instantiations: O(1)
 * * Other operations: O(1)
 *
 * \include test/snippet/utility/type_pack/pack_traits_size.cpp
 */
template <typename... pack_t>
inline constexpr size_t size = sizeof...(pack_t);

/*!\brief Count the occurrences of a type in a pack.
 * \tparam query_t  The type you are searching for.
 * \tparam pack_t   The type pack.
 * \returns The number of occurrences of the `query_t` in `pack_t`.
 * \ingroup utility_type_pack
 *
 * \details
 *
 * ### (Compile-time) Complexity
 *
 * * Number of template instantiations: O(1)
 * * Other operations: O(n)
 *
 * \include test/snippet/utility/type_pack/pack_traits_count.cpp
 */
template <typename query_t, typename... pack_t>
inline constexpr ptrdiff_t count = (std::is_same_v<query_t, pack_t> + ... + 0);

/*!\brief Get the index of the first occurrence of a type in a pack.
 * \tparam query_t  The type you are searching for.
 * \tparam pack_t   The type pack.
 * \returns The position of the first occurrence of `query_t` in `pack_t` or `-1` if it is not contained.
 * \ingroup utility_type_pack
 *
 * \details
 *
 * ### (Compile-time) Complexity
 *
 * * Number of template instantiations: O(1)
 * * Other operations: O(n), possibly == `i`, where `i` is the return value
 *
 * \include test/snippet/utility/type_pack/pack_traits_find.cpp
 */
template <typename query_t, typename... pack_t>
inline constexpr ptrdiff_t find = seqan3::pack_traits::detail::find<query_t, pack_t...>();

/*!\brief Get the index of the first type in a pack that satisfies the given predicate.
 * \tparam pred_t   The predicate that is being evaluated (a class template).
 * \tparam pack_t   The type pack.
 * \returns The index or `-1` if no types match.
 * \ingroup utility_type_pack
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
 * \include test/snippet/utility/type_pack/pack_traits_find_if.cpp
 */
template <template <typename> typename pred_t, typename... pack_t>
inline constexpr ptrdiff_t find_if = seqan3::pack_traits::detail::find_if<pred_t, pack_t...>();

/*!\brief Whether a type occurs in a pack or not.
 * \tparam query_t  The type you are searching for.
 * \tparam pack_t   The type pack.
 * \returns `true` or `false`.
 * \ingroup utility_type_pack
 *
 * \details
 *
 * ### (Compile-time) Complexity
 *
 * * Number of template instantiations: O(1)
 * * Other operations: O(n), possibly == `i`, where `i` is the index of the first occurrence
 *
 * \include test/snippet/utility/type_pack/pack_traits_find.cpp
 */
template <typename query_t, typename... pack_t>
inline constexpr bool contains = (find<query_t, pack_t...> != -1);
//!\}

/*!\name Type pack traits (return a single type)
 * \{
 */

/*!\brief Return the type at given index from the type pack.
 * \tparam idx    The index; must be smaller than the size of the pack.
 * \tparam pack_t The type pack.
 * \ingroup utility_type_pack
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
 * \include test/snippet/utility/type_pack/pack_traits_at.cpp
 */
template <ptrdiff_t idx, typename... pack_t>
    requires (idx >= 0 && idx < sizeof...(pack_t)) || (-idx <= sizeof...(pack_t))
using at = typename decltype(detail::at<idx, pack_t...>())::type;

/*!\brief Return the first type from the type pack.
 * \tparam pack_t The type pack.
 * \ingroup utility_type_pack
 *
 * \details
 *
 * ### (Compile-time) Complexity
 *
 * * Number of template instantiations: O(1)
 * * Other operations: O(1)
 *
 * \include test/snippet/utility/type_pack/pack_traits_front.cpp
 */
template <typename... pack_t>
    requires (sizeof...(pack_t) > 0)
using front = typename decltype(detail::front<pack_t...>())::type;

/*!\brief Return the last type from the type pack.
 * \tparam pack_t The type pack.
 * \ingroup utility_type_pack
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
 * \include test/snippet/utility/type_pack/pack_traits_back.cpp
 */
template <typename... pack_t>
    requires (sizeof...(pack_t) > 0)
using back = typename decltype((std::type_identity<pack_t>{}, ...))::type; // use comma operator

//!\}

/*!\name Type pack traits (return a type list)
 * \{
 */

/*!\brief Return a seqan3::type_list of all the types in the type pack, except the first.
 * \tparam pack_t The type pack.
 * \ingroup utility_type_pack
 *
 * \details
 *
 * ### (Compile-time) Complexity
 *
 * * Number of template instantiations: O(1)
 * * Other operations: O(1)
 *
 * \include test/snippet/utility/type_pack/pack_traits_drop_front.cpp
 */
template <typename... pack_t>
    requires (sizeof...(pack_t) > 0)
using drop_front = typename decltype(detail::drop_front<pack_t...>())::type;

/*!\brief Apply a transformation trait to every type in the pack and return a seqan3::type_list of the results.
 * \tparam trait_t The transformation trait, **can be an alias template**, e.g. a transformation trait shortcut.
 * \tparam pack_t  The type pack.
 * \ingroup utility_type_pack
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
 * \include test/snippet/utility/type_pack/pack_traits_transform.cpp
 */
template <template <typename> typename trait_t, typename... pack_t>
using transform = seqan3::type_list<trait_t<pack_t>...>;

//!\}

/*!\name Type pack traits (return a type list)
 * \{
 */

/*!\brief Return a seqan3::type_list of the first `n` types in the type pack.
 * \tparam i        The target size; must be >= 0 and <= the size of the type pack.
 * \tparam pack_t   The type pack.
 * \ingroup utility_type_pack
 *
 * \details
 *
 * ### (Compile-time) Complexity
 *
 * * Number of template instantiations: O(n)
 * * Other operations: O(n)
 *
 * \include test/snippet/utility/type_pack/pack_traits_take.cpp
 */
template <ptrdiff_t i, typename... pack_t>
    requires (i >= 0 && i <= size<pack_t...>)
using take = typename decltype(detail::split_after<i, pack_t...>(type_list<>{}))::first_type;

/*!\brief Return a seqan3::type_list of the types in the type pack, except the first `n`.
 * \tparam i        The amount to drop; must be >= 0 and <= the size of the type pack.
 * \tparam pack_t   The type pack.
 * \ingroup utility_type_pack
 *
 * \details
 *
 * ### (Compile-time) Complexity
 *
 * * Number of template instantiations: O(n)
 * * Other operations: O(n)
 *
 * \include test/snippet/utility/type_pack/pack_traits_drop.cpp
 */
template <ptrdiff_t i, typename... pack_t>
    requires (i >= 0 && i <= size<pack_t...>)
using drop = typename decltype(detail::split_after<i, pack_t...>(type_list<>{}))::second_type;

/*!\brief Return a seqan3::type_list of the last `n` types in the type pack.
 * \tparam i        The target size; must be >= 0 and <= the size of the type pack.
 * \tparam pack_t   The type pack.
 * \ingroup utility_type_pack
 *
 * \details
 *
 * ### (Compile-time) Complexity
 *
 * * Number of template instantiations: O(n)
 * * Other operations: O(n)
 *
 * \include test/snippet/utility/type_pack/pack_traits_take_last.cpp
 */
template <ptrdiff_t i, typename... pack_t>
    requires (i >= 0 && i <= size<pack_t...>)
using take_last = drop<size<pack_t...> - i, pack_t...>;

/*!\brief Return a seqan3::type_list of the types the type pack, except the last `n`.
 * \tparam i        The amount to drop; must be >= 0 and <= the size of the type pack.
 * \tparam pack_t   The type pack.
 * \ingroup utility_type_pack
 *
 * \details
 *
 * ### (Compile-time) Complexity
 *
 * * Number of template instantiations: O(n)
 * * Other operations: O(n)
 *
 * \include test/snippet/utility/type_pack/pack_traits_take_last.cpp
 */
template <ptrdiff_t i, typename... pack_t>
    requires (i >= 0 && i <= size<pack_t...>)
using drop_last = take<size<pack_t...> - i, pack_t...>;

/*!\brief Split a type pack into two parts returned as a pair of seqan3::type_list.
 * \tparam i        The number of elements after which to split; must be >= 0 and <= the size of the type pack.
 * \tparam pack_t   The type pack.
 * \ingroup utility_type_pack
 *
 * \details
 *
 * ### (Compile-time) Complexity
 *
 * * Number of template instantiations: O(n)
 * * Other operations: O(n)
 *
 * \include test/snippet/utility/type_pack/pack_traits_take_last.cpp
 */
template <ptrdiff_t i, typename... pack_t>
    requires (i >= 0 && i <= size<pack_t...>)
using split_after = decltype(detail::split_after<i, pack_t...>(type_list<>{}));

/*!\brief Replace the type at the given index with the given type.
 * \tparam replace_t The type to replace the old type with.
 * \tparam i         The index of the type to be replaced.
 * \tparam pack_t    The (input) type pack.
 * \ingroup utility_type_pack
 *
 * \details
 *
 * ### (Compile-time) Complexity
 *
 * * Number of template instantiations: O(n)
 * * Other operations: O(n)
 *
 * \include test/snippet/utility/type_pack/pack_traits_take_last.cpp
 */
template <typename replace_t, std::ptrdiff_t i, typename... pack_t>
    requires (i >= 0 && i < size<pack_t...>)
using replace_at = decltype(detail::replace_at<replace_t, i, pack_t...>(std::make_index_sequence<size<pack_t...>>{}));

//!\}

} // namespace seqan3::pack_traits
