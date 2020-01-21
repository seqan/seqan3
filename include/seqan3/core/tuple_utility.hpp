// -----------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2020, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2020, Knut Reinert & MPI für molekulare Genetik
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at: https://github.com/seqan/seqan3/blob/master/LICENSE.md
// -----------------------------------------------------------------------------------------------------

/*!\file
 * \brief Provides utility functions for tuple like interfaces.
 * \author Rene Rahn <rene.rahn AT fu-berlin.de>
 */

#pragma once

#include <utility>

#include <seqan3/core/concept/tuple.hpp>
#include <seqan3/core/pod_tuple.hpp>
#include <seqan3/core/type_list/traits.hpp>

namespace seqan3::detail
{

/*!\brief Helper function for seqan3::tuple_split.
 * \ingroup core
 *
 * \tparam    beg     A template value containing the start position from where to extract the values.
 * \tparam    tuple_t A template alias for a tuple like object.
 * \tparam    ...ts   Types `tuple_t` is specified with.
 * \tparam    ...Is   Indices of the tuple elements that should be extracted.
 *
 * \param[in] t       The original tuple to split.
 * \param[in] idx     An std::index_sequence with all indices that should be extracted beginning at `beg`.
 *
 * \returns A new tuple with the extracted elements.
 */
template <size_t beg,
          template <typename ...> typename tuple_t,
          size_t ... Is,
          typename ...ts>
//!\cond
    requires tuple_like<tuple_t<ts...>> && tuple_like<tuple_t<>>
//!\endcond
constexpr auto tuple_split(tuple_t<ts...> const & t, std::index_sequence<Is...> const & SEQAN3_DOXYGEN_ONLY(idx))
{
    return tuple_t<std::tuple_element_t<beg + Is, tuple_t<ts...>>...>{std::get<beg + Is>(t)...};
}

//!\copydoc seqan3::detail::tuple_split
template <size_t beg,
          template <typename ...> typename tuple_t,
          size_t ... Is,
          typename ...ts>
//!\cond
    requires tuple_like<tuple_t<ts...>> && tuple_like<tuple_t<>>
//!\endcond
constexpr auto tuple_split(tuple_t<ts...> && t, std::index_sequence<Is...> const & SEQAN3_DOXYGEN_ONLY(idx))
{
    return tuple_t<std::tuple_element_t<beg + Is, tuple_t<ts...>>...>{std::move(std::get<beg + Is>(t))...};
}
} // namespace seqan3::detail

namespace seqan3
{
/*!\name Tuple utility functions
 * \brief Helper functions for tuple like objects.
 * \{
 */
/*!\brief Splits a tuple like data structure at the given position.
 * \ingroup core
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
 * \include test/snippet/core/tuple_utility.cpp
 *
 * ### Complexity
 *
 * Linear in the number of elements.
 *
 * ### Thread safety
 *
 * Concurrent invocations of this functions are thread safe.
 */
template <size_t pivot_c, template <typename ...> typename tuple_t, typename ...ts>
//!\cond
    requires tuple_like<tuple_t<ts...>>
//!\endcond
constexpr auto tuple_split(tuple_t<ts...> const & t)
{
    static_assert(pivot_c <= sizeof...(ts));

    return tuple_t{detail::tuple_split<0>(t, std::make_index_sequence<pivot_c>{}),
                   detail::tuple_split<pivot_c>(t, std::make_index_sequence<sizeof...(ts) - pivot_c>{})};
}

//!\copydoc seqan3::tuple_split
template <size_t pivot_c, template <typename ...> typename tuple_t, typename ...ts>
//!\cond
    requires tuple_like<tuple_t<ts...>>
//!\endcond
constexpr auto tuple_split(tuple_t<ts...> && t)
{
    static_assert(pivot_c <= sizeof...(ts));

    return tuple_t{detail::tuple_split<0>(std::move(t), std::make_index_sequence<pivot_c>{}),
                   detail::tuple_split<pivot_c>(std::move(t), std::make_index_sequence<sizeof...(ts) - pivot_c>{})};
}

/*!\brief Splits a tuple like data structure at the first position of the given type.
 * \ingroup core
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
 * \include test/snippet/core/tuple_utility.cpp
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
    constexpr size_t pivot_c = list_traits::find<pivot_t, detail::tuple_type_list_t<remove_cvref_t<tuple_t>>>;

    static_assert(pivot_c <= std::tuple_size_v<remove_cvref_t<tuple_t>>);

    return tuple_split<pivot_c>(std::forward<tuple_t>(t));
}

/*!\brief Removes the first element of a tuple.
 * \ingroup core
 *
 * \param[in] t  The original tuple.
 *
 * \returns A new tuple without the first element of `t`.
 *
 * \details
 *
 * Note, that the tuple must contain at least one element and must support empty tuple types, i.e. std::pair cannot
 * be used.
 *
 * ### Complexity
 *
 * Linear in the number of elements.
 *
 * ### Thread safety
 *
 * Concurrent invocations of this functions are thread safe.
 */
template <tuple_like tuple_t>
constexpr auto tuple_pop_front(tuple_t && t)
{
    static_assert(std::tuple_size_v<remove_cvref_t<tuple_t>> > 0);

    return std::get<1>(tuple_split<1>(std::forward<tuple_t>(t)));
}
//!\}

} // namespace seqan3
