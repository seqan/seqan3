// SPDX-FileCopyrightText: 2006-2025 Knut Reinert & Freie Universität Berlin
// SPDX-FileCopyrightText: 2016-2025 Knut Reinert & MPI für molekulare Genetik
// SPDX-License-Identifier: BSD-3-Clause

/*!\file
 * \brief Provides seqan3::tuple_like.
 * \author Hannes Hauswedell <hannes.hauswedell AT fu-berlin.de>
 */

#pragma once

#include <concepts>
#include <tuple>
#include <type_traits>

#include <seqan3/core/detail/template_inspection.hpp>
#include <seqan3/utility/tuple/pod_tuple.hpp>
#include <seqan3/utility/type_list/type_list.hpp>
#include <seqan3/utility/type_traits/basic.hpp>

namespace seqan3::detail
{

/*!\interface seqan3::detail::tuple_size <>
 * \ingroup utility_tuple
 * \brief Subconcept definition for seqan3::tuple_like to test for std::tuple_size-interface.
 * \see seqan3::tuple_like
 */
//!\cond
template <typename tuple_t>
concept tuple_size = requires (tuple_t v) {
    { std::tuple_size<tuple_t>::value } -> std::convertible_to<size_t>;
};
//!\endcond

/*!\interface seqan3::detail::tuple_get <>
 * \ingroup utility_tuple
 * \brief Subconcept definition for seqan3::tuple_like to test for std::get-interface.
 * \see seqan3::tuple_like
 */
//!\cond
template <typename tuple_t>
concept tuple_get = requires (tuple_t & v, tuple_t const & v_c) {
    requires std::tuple_size_v<tuple_t> > 0;

    typename std::tuple_element<0, tuple_t>::type;

    { get<0>(v) } -> std::convertible_to<typename std::tuple_element<0, tuple_t>::type>;
    //     requires weakly_assignable_from<decltype(get<0>(v)), typename std::tuple_element<0, tuple_t>::type>;
    //TODO check that the previous returns something that can be assigned to
    // unfortunately std::assignable_from requires lvalue-reference, but we want to accept xvalues too (returned
    // proxies)
    { get<0>(v_c) } -> std::convertible_to<typename std::tuple_element<0, tuple_t>::type>;
    { get<0>(std::move(v)) } -> std::convertible_to<typename std::tuple_element<0, tuple_t>::type>;
    { get<0>(std::move(v_c)) } -> std::convertible_to<typename std::tuple_element<0, tuple_t>::type const &&>;
};
//!\endcond

/*!\brief Transformation trait to expose the tuple element types as seqan3::type_list
 * \ingroup utility_tuple
 * \tparam tuple_t The tuple to extract the element types from.
 *
 * \returns A seqan3::type_list over the element types of the given tuple.
 * \see seqan3::detail::tuple_type_list_t
 */
template <detail::tuple_size tuple_t>
struct tuple_type_list
{
protected:
    //!\brief Helper function to extract the types using the tuple elements.
    template <size_t... Is>
    static constexpr auto invoke_to_type_list(std::index_sequence<Is...>)
    {
        return type_list<std::tuple_element_t<Is, tuple_t>...>{};
    }

public:
    //!\brief The generated seqan3::type_list.
    using type = decltype(invoke_to_type_list(std::make_index_sequence<std::tuple_size<tuple_t>::value>{}));
};

/*!\brief Helper type for seqan3::detail::tuple_type_list
 * \ingroup utility_tuple
 *
 * \see seqan3::detail::tuple_type_list
 */
template <detail::tuple_size tuple_t>
using tuple_type_list_t = typename tuple_type_list<tuple_t>::type;

/*!\brief Helper type function to check for std::totally_ordered on all elements of the given tuple type.
 * \ingroup utility_tuple
 */
template <typename... elements_t>
inline constexpr auto all_elements_model_totally_ordered(seqan3::type_list<elements_t...>)
    -> std::bool_constant<(std::totally_ordered<elements_t> && ... && true)>;

/*!\brief Helper type trait function to check for std::totally_ordered on all elements of the given tuple type.
 * \tparam tuple_t The tuple to check if all elements model std::totally_ordered.
 * \ingroup utility_tuple
 */
template <typename tuple_t>
    requires requires () {
        { detail::all_elements_model_totally_ordered(tuple_type_list_t<tuple_t>{}) };
    }
static constexpr bool all_elements_model_totally_ordered_v =
    decltype(detail::all_elements_model_totally_ordered(tuple_type_list_t<tuple_t>{}))::value;
} // namespace seqan3::detail

namespace seqan3
{

// ----------------------------------------------------------------------------
// tuple_like
// ----------------------------------------------------------------------------

/*!\interface seqan3::tuple_like
 * \extends std::totally_ordered
 * \ingroup utility_tuple
 * \brief Whether a type behaves like a tuple.
 *
 * \details
 *
 * Types that meet this concept are for example std::tuple, std::pair, std::array, seqan3::pod_tuple, seqan3::record.
 * The std::totally_ordered will only be required if all types contained in the tuple-like
 * data structure are themselves strict totally ordered.
 */
/*!\name Requirements for seqan3::tuple_like
 * \brief You can expect these (meta-)functions on all types that implement seqan3::tuple_like.
 * \{
 */
/*!\var         size_t std::tuple_size_v<type>
 * \brief       A unary type trait that holds the number of elements in the tuple.
 * \tparam      type The tuple-like type.
 * \relates     seqan3::tuple_like
 *
 * \details
 * \attention This is a concept requirement, not an actual function (however types satisfying this concept
 * will provide an implementation).
 */
/*!\typedef     std::tuple_elment_t<i, type>
 * \brief       A transformation trait that holds the type of elements in the tuple.
 * \tparam      i Index of the queried element type.
 * \tparam      type The tuple-like type.
 * \relates     seqan3::tuple_like
 *
 * \details
 * \attention This is a concept requirement, not an actual function (however types satisfying this concept
 * will provide an implementation).
 * \attention This constraint is **not enforced** since empty tuples are valid.
 */
/*!\fn              auto && std::get<i>(type && val)
 * \brief           Return the i-th element of the tuple.
 * \relates         seqan3::tuple_like
 * \tparam          i The index of the element to return (of type `size_t`).
 * \param[in,out]   val The tuple-like object to operate on.
 * \returns         The i-th value in the tuple.
 *
 * \details
 * \attention This is a concept requirement, not an actual function (however types satisfying this concept
 * will provide an implementation).
 * \attention This constraint is **not enforced** since empty tuples are valid.
 */
//!\}
//!\cond
template <typename t>
concept tuple_like = detail::tuple_size<std::remove_reference_t<t>> && requires (t v) {
    typename detail::tuple_type_list<std::remove_cvref_t<t>>::type;

    // NOTE(rrahn): To check the full tuple_concept including the get interface and the std::totally_ordered
    //              we need to make some assumptions. In general these checks can only be executed if the tuple is not
    //              empty. Furthermore, the std::totally_ordered can only be checked if all elements in the
    //              tuple are strict_totally_ordered. This is done, by the fold expression in the second part.
    requires (std::tuple_size<std::remove_reference_t<t>>::value == 0)
                 || (detail::tuple_get<std::remove_cvref_t<t>>
                     && (!detail::all_elements_model_totally_ordered_v<std::remove_cvref_t<t>>
                         || std::totally_ordered<std::remove_cvref_t<t>>));
};
//!\endcond

/*!\interface seqan3::pair_like
 * \extends seqan3::tuple_like
 * \ingroup utility_tuple
 * \brief Whether a type behaves like a tuple with exactly two elements.
 *
 * \details
 *
 * Types that meet this concept are for example std::tuple, std::pair, std::array, seqan3::pod_tuple,
 * iff std::tuple_size equals `2`.
 */
//!\cond
template <typename t>
concept pair_like = tuple_like<t> && std::tuple_size_v<std::remove_reference_t<t>> == 2;
//!\endcond

} // namespace seqan3
