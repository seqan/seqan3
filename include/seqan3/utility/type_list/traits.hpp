// SPDX-FileCopyrightText: 2006-2025 Knut Reinert & Freie Universität Berlin
// SPDX-FileCopyrightText: 2016-2025 Knut Reinert & MPI für molekulare Genetik
// SPDX-License-Identifier: BSD-3-Clause

/*!\file
 * \author Hannes Hauswedell <hannes.hauswedell AT fu-berlin.de>
 * \brief Provides traits for seqan3::type_list.
 */

#pragma once

#include <type_traits>

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
 * \ingroup utility_type_list
 */
template <ptrdiff_t idx, typename... pack_t>
std::type_identity<seqan3::pack_traits::at<idx, pack_t...>> at(type_list<pack_t...>);

/*!\brief Implementation for seqan3::list_traits::front.
 * \tparam pack_t   Types in the type list.
 * \ingroup utility_type_list
 */
template <typename... pack_t>
std::type_identity<seqan3::pack_traits::front<pack_t...>> front(type_list<pack_t...>);

/*!\brief Implementation for seqan3::list_traits::back.
 * \tparam pack_t   Types in the type list.
 * \ingroup utility_type_list
 */
template <typename... pack_t>
std::type_identity<seqan3::pack_traits::back<pack_t...>> back(type_list<pack_t...>);

/*!\brief Implementation for seqan3::list_traits::concat.
 * \tparam pack1_t   Types in the first type list.
 * \tparam pack2_t   Types in the second type list.
 * \ingroup utility_type_list
 */
template <typename... pack1_t, typename... pack2_t>
type_list<pack1_t..., pack2_t...> concat(type_list<pack1_t...>, type_list<pack2_t...>);

/*!\brief Implementation for seqan3::list_traits::concat [overload for more than two lists].
 * \tparam pack1_t      Types in the first type list.
 * \tparam pack2_t      Types in the second type list.
 * \tparam more_lists_t The remaining type lists.
 * \ingroup utility_type_list
 */
template <typename... pack1_t, typename... pack2_t, typename... more_lists_t>
auto concat(type_list<pack1_t...>, type_list<pack2_t...>, more_lists_t...)
{
    return concat(type_list<pack1_t..., pack2_t...>{}, more_lists_t{}...);
}

/*!\brief Implementation for seqan3::list_traits::drop_front.
 * \tparam pack_t   Types in the type list.
 * \ingroup utility_type_list
 */
template <typename... pack_t>
pack_traits::drop_front<pack_t...> drop_front(type_list<pack_t...>);

/*!\brief Implementation for seqan3::list_traits::transform.
 * \tparam trait_t The trait to transform, **must be an alias template**, e.g. a transformation trait shortcut.
 * \tparam pack_t  Types in the type list.
 * \ingroup utility_type_list
 */
template <template <typename> typename trait_t, typename... pack_t>
pack_traits::transform<trait_t, pack_t...> transform(type_list<pack_t...>);

/*!\brief Implementation for seqan3::list_traits::split_after.
 * \tparam idx The index after which to split.
 * \tparam pack_t Types in the type list to split
 * \ingroup utility_type_list
 */
template <ptrdiff_t idx, typename... pack1_t>
pack_traits::split_after<idx, pack1_t...> split_after(type_list<pack1_t...>);

/*!\brief Implementation for seqan3::list_traits::repeat.
 * \tparam count The number of repititions.
 * \tparam t The type to repeat
 * \ingroup utility_type_list
 */
template <size_t count, typename t>
auto repeat()
{
    if constexpr (count == 0)
        return type_list<>{};
    else if constexpr (count == 1)
        return type_list<t>{};
    else if constexpr (count == 2)
        return type_list<t, t>{};
    else if constexpr (count == 3)
        return type_list<t, t, t>{};
    else if constexpr (count == 4)
        return type_list<t, t, t, t>{};
    else if constexpr (count == 5)
        return type_list<t, t, t, t, t>{};
    else
        return concat(repeat<5, t>(), repeat<count - 5, t>());
}

/*!\brief Implementation for seqan3::list_traits::replace_at.
 * \tparam replace_t The type replacing the old one.
 * \tparam idx The index of the type to replace.
 * \tparam pack_t Types in the type list to be modified.
 * \ingroup utility_type_list
 */
template <typename replace_t, ptrdiff_t idx, typename... pack_t>
pack_traits::replace_at<replace_t, idx, pack_t...> replace_at(type_list<pack_t...>);

//!\brief A replacement for meta::reverse [recursion anchor]
inline constexpr type_list<> reverse(type_list<>)
{
    return {};
}

//!\brief A replacement for meta::reverse [recursion]
template <typename head_t, typename... pack_t>
auto reverse(type_list<head_t, pack_t...>)
{
    return concat(reverse(type_list<pack_t...>{}), type_list<head_t>{});
}

//!\brief Constructs the multiset difference `list1 \ list2` [recursion anchor]
template <typename... current_list_t>
constexpr seqan3::type_list<current_list_t...> type_list_difference(seqan3::type_list<current_list_t...>,
                                                                    seqan3::type_list<>)
{
    return {};
}

//!\brief Constructs the multiset difference `list1 \ list2` [recursion]
template <typename... current_list_t, typename remove_t, typename... remove_list_t>
constexpr auto type_list_difference(seqan3::type_list<current_list_t...>, seqan3::type_list<remove_t, remove_list_t...>)
{
    constexpr auto pos = seqan3::pack_traits::find<remove_t, current_list_t...>;
    if constexpr (pos >= 0)
    {
        using split_list_t = seqan3::pack_traits::split_after<pos, current_list_t...>;

        using split_list1_t = typename split_list_t::first_type;
        using split_list2_t = decltype(drop_front(typename split_list_t::second_type{}));
        using filtered_list_t = decltype(concat(split_list1_t{}, split_list2_t{}));
        return type_list_difference(filtered_list_t{}, seqan3::type_list<remove_t, remove_list_t...>{});
    }
    else
    {
        // remove_t not contained in current_list_t
        using filtered_list_t = seqan3::type_list<current_list_t...>;
        return type_list_difference(filtered_list_t{}, seqan3::type_list<remove_list_t...>{});
    }
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
 * \ingroup utility_type_list
 * \copydetails seqan3::pack_traits::size
 *
 * \include test/snippet/utility/type_list/list_traits_size.cpp
 */
template <typename... pack_t>
inline constexpr size_t size<type_list<pack_t...>> = sizeof...(pack_t);

//!\cond
template <typename query_t, typename list_t>
    requires seqan3::detail::template_specialisation_of<list_t, seqan3::type_list>
inline constexpr ptrdiff_t count = -1;
//!\endcond

/*!\brief Count the occurrences of a type in a type list.
 * \ingroup utility_type_list
 * \copydetails seqan3::pack_traits::count
 *
 * \include test/snippet/utility/type_list/list_traits_count.cpp
 */
template <typename query_t, typename... pack_t>
inline constexpr ptrdiff_t count<query_t, type_list<pack_t...>> = seqan3::pack_traits::count<query_t, pack_t...>;

//!\cond
template <typename query_t, typename list_t>
    requires seqan3::detail::template_specialisation_of<list_t, seqan3::type_list>
inline constexpr ptrdiff_t find = -1;
//!\endcond

/*!\brief Get the index of the first occurrence of a type in a type list.
 * \ingroup utility_type_list
 * \copydetails seqan3::pack_traits::find
 *
 * \include test/snippet/utility/type_list/list_traits_find.cpp
 */
template <typename query_t, typename... pack_t>
inline constexpr ptrdiff_t find<query_t, type_list<pack_t...>> =
    seqan3::pack_traits::detail::find<query_t, pack_t...>();

//!\cond
template <template <typename> typename pred_t, typename list_t>
    requires seqan3::detail::template_specialisation_of<list_t, seqan3::type_list>
inline constexpr ptrdiff_t find_if = -1;
//!\endcond

/*!\brief Get the index of the first type in a type list that satisfies the given predicate.
 * \ingroup utility_type_list
 * \copydetails seqan3::pack_traits::find_if
 *
 * \include test/snippet/utility/type_list/list_traits_find_if.cpp
 */
template <template <typename> typename pred_t, typename... pack_t>
inline constexpr ptrdiff_t find_if<pred_t, type_list<pack_t...>> =
    seqan3::pack_traits::detail::find_if<pred_t, pack_t...>();

/*!\brief Whether a type occurs in a type list or not.
 * \ingroup utility_type_list
 * \copydetails seqan3::pack_traits::contains
 *
 * \include test/snippet/utility/type_list/list_traits_contains.cpp
 */
template <typename query_t, typename list_t>
    requires seqan3::detail::template_specialisation_of<list_t, seqan3::type_list>
inline constexpr bool contains = (find<query_t, list_t> != -1);

//!\}

/*!\name Type list traits (return a single type)
 * \{
 */

/*!\brief Return the type at given index from the type list.
 * \tparam idx    The index; must be smaller than the size of the type list.
 * \tparam list_t The type_list.
 * \ingroup utility_type_list
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
    requires (seqan3::detail::template_specialisation_of<list_t, seqan3::type_list>)
              && ((idx >= 0 && idx < size<list_t>) || (-idx <= size<list_t>))
using at = typename decltype(detail::at<idx>(list_t{}))::type;

/*!\brief Return the first type from the type list.
 * \tparam list_t The type list.
 * \ingroup utility_type_list
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
    requires (seqan3::detail::template_specialisation_of<list_t, seqan3::type_list>) && (size<list_t> > 0)
using front = typename decltype(detail::front(list_t{}))::type;

/*!\brief Return the last type from the type list.
 * \tparam list_t The type list.
 * \ingroup utility_type_list
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
    requires (seqan3::detail::template_specialisation_of<list_t, seqan3::type_list>) && (size<list_t> > 0)
using back = typename decltype(detail::back(list_t{}))::type;

//!\}

/*!\name Type list traits (return a type list)
 * \{
 */

/*!\brief Join two seqan3::type_list s into one.
 * \tparam list1_t The first (input) type list.
 * \tparam list2_t The second (input) type list.
 * \ingroup utility_type_list
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
template <typename... lists_t>
    requires (seqan3::detail::template_specialisation_of<lists_t, seqan3::type_list> && ...)
using concat = decltype(detail::concat(lists_t{}...));

/*!\brief Return a seqan3::type_list of all the types in the type list, except the first.
 * \tparam list_t The (input) type list.
 * \ingroup utility_type_list
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
    requires (seqan3::detail::template_specialisation_of<list_t, seqan3::type_list>) && (size<list_t> > 0)
using drop_front = decltype(detail::drop_front(list_t{}));

/*!\brief Return a seqan3::type_list of the first `n` types in the input type list.
 * \tparam i        The target size; must be >= 0 and <= the size of the input type list.
 * \tparam list_t   The (input) type list.
 * \ingroup utility_type_list
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
    requires (seqan3::detail::template_specialisation_of<list_t, seqan3::type_list>) && (i >= 0 && i <= size<list_t>)
using take = typename decltype(detail::split_after<i>(list_t{}))::first_type;

/*!\brief Return a seqan3::type_list of the types in the input type list, except the first `n`.
 * \tparam i        The amount to drop; must be >= 0 and <= the size of the input type list.
 * \tparam list_t   The (input) type list.
 * \ingroup utility_type_list
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
    requires (seqan3::detail::template_specialisation_of<list_t, seqan3::type_list>) && (i >= 0 && i <= size<list_t>)
using drop = typename decltype(detail::split_after<i>(list_t{}))::second_type;

/*!\brief Return a seqan3::type_list of the last `n` types in the input type list.
 * \tparam i        The target size; must be >= 0 and <= the size of the input type list.
 * \tparam list_t   The (input) type list.
 * \ingroup utility_type_list
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
    requires (seqan3::detail::template_specialisation_of<list_t, seqan3::type_list>) && (i >= 0 && i <= size<list_t>)
using take_last = drop<size<list_t> - i, list_t>;

/*!\brief Return a seqan3::type_list of the types the input type list, except the last `n`.
 * \tparam i        The amount to drop; must be >= 0 and <= the size of the input type list.
 * \tparam list_t   The (input) type list.
 * \ingroup utility_type_list
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
    requires (seqan3::detail::template_specialisation_of<list_t, seqan3::type_list>) && (i >= 0 && i <= size<list_t>)
using drop_last = take<size<list_t> - i, list_t>;

/*!\brief Split a seqan3::type_list into two parts returned as a pair of seqan3::type_list.
 * \tparam i        The number of elements after which to split; must be >= 0 and <= the size of the input type list.
 * \tparam list_t   The (input) type list.
 * \ingroup utility_type_list
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
    requires (seqan3::detail::template_specialisation_of<list_t, seqan3::type_list>) && (i >= 0 && i <= size<list_t>)
using split_after = decltype(detail::split_after<i>(list_t{}));

/*!\brief Apply a transformation trait to every type in the list and return a seqan3::type_list of the results.
 * \tparam trait_t The trait to transform, **must be an alias template**, e.g. a transformation trait shortcut.
 * \tparam list_t  The (input) type list.
 * \ingroup utility_type_list
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
    requires seqan3::detail::template_specialisation_of<list_t, seqan3::type_list>
using transform = decltype(detail::transform<trait_t>(list_t{}));

/*!\brief Replace the type at the given index with the given type.
 * \tparam replace_t The type to replace the old type with.
 * \tparam i         The index of the type to be replaced.
 * \tparam list_t    The (input) type list.
 * \ingroup utility_type_list
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
    requires (seqan3::detail::template_specialisation_of<list_t, seqan3::type_list>) && (i >= 0 && i < size<list_t>)
using replace_at = decltype(detail::replace_at<replace_t, i>(list_t{}));

/*!\brief Create a type list with the given type repeated `count` times..
 * \tparam count The number of repititions.
 * \tparam t The type to repeat
 * \ingroup utility_type_list
 *
 * ### (Compile-time) Complexity
 *
 * * Number of template instantiations: O(n)
 * * Other operations: O(n)
 *
 * However, with a constant of 0.2, e.g. repeat<10,int> results in 2-3 instantiations.
 */
template <size_t count, typename t>
using repeat = decltype(detail::repeat<count, t>());
//!\}

} // namespace seqan3::list_traits
