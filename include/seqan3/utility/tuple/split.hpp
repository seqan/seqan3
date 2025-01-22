// SPDX-FileCopyrightText: 2006-2025 Knut Reinert & Freie Universität Berlin
// SPDX-FileCopyrightText: 2016-2025 Knut Reinert & MPI für molekulare Genetik
// SPDX-License-Identifier: BSD-3-Clause

/*!\file
 * \brief Provides seqan3::tuple_split.
 * \author Rene Rahn <rene.rahn AT fu-berlin.de>
 */

#pragma once

#include <utility>

#include <seqan3/utility/tuple/concept.hpp>
#include <seqan3/utility/type_list/traits.hpp>

namespace seqan3::detail
{

/*!\brief Helper function for seqan3::tuple_split.
 * \ingroup utility_tuple
 *
 * \tparam    beg     A template value containing the start position from where to extract the values.
 * \tparam    tuple_t A template alias for a tuple like object.
 * \tparam    ...ts   Types `tuple_t` is specified with.
 * \tparam    ...Is   Indices of the tuple elements that should be extracted.
 *
 * \param[in] t       The original tuple to split.
 * \param[in] idx     A std::index_sequence with all indices that should be extracted beginning at `beg`.
 *
 * \returns A new tuple with the extracted elements.
 */
template <size_t beg, template <typename...> typename tuple_t, size_t... Is, typename... ts>
    requires tuple_like<tuple_t<ts...>> && tuple_like<tuple_t<>>
constexpr auto tuple_split(tuple_t<ts...> const & t, std::index_sequence<Is...> const & SEQAN3_DOXYGEN_ONLY(idx))
{
#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wmissing-braces"
    return tuple_t<std::tuple_element_t<beg + Is, tuple_t<ts...>>...>{std::get<beg + Is>(t)...};
#pragma GCC diagnostic pop
}

//!\copydoc seqan3::detail::tuple_split
template <size_t beg, template <typename...> typename tuple_t, size_t... Is, typename... ts>
    requires tuple_like<tuple_t<ts...>> && tuple_like<tuple_t<>>
constexpr auto tuple_split(tuple_t<ts...> && t, std::index_sequence<Is...> const & SEQAN3_DOXYGEN_ONLY(idx))
{
#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wmissing-braces"
    return tuple_t<std::tuple_element_t<beg + Is, tuple_t<ts...>>...>{std::move(std::get<beg + Is>(t))...};
#pragma GCC diagnostic pop
}
} // namespace seqan3::detail

namespace seqan3
{
/*!\brief Splits a tuple like data structure at the given position.
 * \ingroup utility_tuple
 *
 * \tparam    pivot_c A template value specifying the split position.
 * \tparam    tuple_t A template alias for a tuple like object.
 * \tparam    ...ts   Types `tuple_t` is specified with.
 * \param[in] t       The original tuple to split.
 *
 * \returns A new tuple of tuples with the left side of the split and the right side of the split.
 *
 * \details
 *
 * Splits a tuple into two tuples, while the element at the split position will be contained in the second tuple.
 * Note, that the returned tuples can be empty. For this reason it is not possible to use tuple like objects,
 * that cannot be empty, i.e. std::pair. Using such an object will emit an compiler error.
 *
 * ### example
 *
 * \include test/snippet/utility/tuple_utility.cpp
 *
 * ### Complexity
 *
 * Linear in the number of elements.
 *
 * ### Thread safety
 *
 * Concurrent invocations of this functions are thread safe.
 */
template <size_t pivot_c, template <typename...> typename tuple_t, typename... ts>
    requires tuple_like<tuple_t<ts...>>
constexpr auto tuple_split(tuple_t<ts...> const & t)
{
    static_assert(pivot_c <= sizeof...(ts));

#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wmissing-braces"
    return tuple_t{detail::tuple_split<0>(t, std::make_index_sequence<pivot_c>{}),
                   detail::tuple_split<pivot_c>(t, std::make_index_sequence<sizeof...(ts) - pivot_c>{})};
#pragma GCC diagnostic pop
}

//!\copydoc seqan3::tuple_split
template <size_t pivot_c, template <typename...> typename tuple_t, typename... ts>
    requires tuple_like<tuple_t<ts...>>
constexpr auto tuple_split(tuple_t<ts...> && t)
{
    static_assert(pivot_c <= sizeof...(ts));
#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wmissing-braces"
    return tuple_t{detail::tuple_split<0>(std::move(t), std::make_index_sequence<pivot_c>{}),
                   detail::tuple_split<pivot_c>(std::move(t), std::make_index_sequence<sizeof...(ts) - pivot_c>{})};
#pragma GCC diagnostic pop
}

/*!\brief Splits a tuple like data structure at the first position of the given type.
 * \ingroup utility_tuple
 *
 * \tparam    pivot_t A template type specifying the split position.
 * \param[in] t       The original tuple to split.
 *
 * \returns A new tuple of tuples with the left side of the split and the right side of the split.
 *
 * \details
 *
 * Splits a tuple into two tuples, while the element at the split position will be contained in the second tuple.
 * Note, that the returned tuples can be empty. For this reason it is not possible to use tuple like objects,
 * that cannot be empty, i.e. std::pair. Using such an object will emit an compiler error.
 *
 * ### example
 *
 * \include test/snippet/utility/tuple_utility.cpp
 *
 * ### Complexity
 *
 * Linear in the number of elements.
 *
 * ### Thread safety
 *
 * Concurrent invocations of this functions are thread safe.
 */
template <typename pivot_t, tuple_like tuple_t>
constexpr auto tuple_split(tuple_t && t)
{
    constexpr size_t pivot_c = list_traits::find<pivot_t, detail::tuple_type_list_t<std::remove_cvref_t<tuple_t>>>;

    static_assert(pivot_c <= std::tuple_size_v<std::remove_cvref_t<tuple_t>>);

    return tuple_split<pivot_c>(std::forward<tuple_t>(t));
}

} // namespace seqan3
