// -----------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2020, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2020, Knut Reinert & MPI für molekulare Genetik
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at: https://github.com/seqan/seqan3/blob/master/LICENSE.md
// -----------------------------------------------------------------------------------------------------

/*!\file
 * \brief Provides seqan3::tuple_like.
 * \author Hannes Hauswedell <hannes.hauswedell AT fu-berlin.de>
 */

#pragma once

#include <tuple>
#include <type_traits>

#include <seqan3/core/type_traits/basic.hpp>
#include <seqan3/core/type_traits/template_inspection.hpp>
#include <seqan3/core/type_list/type_list.hpp>
#include <seqan3/core/pod_tuple.hpp>
#include <seqan3/std/concepts>

namespace seqan3::detail
{

/*!\interface seqan3::detail::tuple_size <>
 * \ingroup   core
 * \brief     Subconcept definition for seqan3::tuple_like to test for std::tuple_size-interface.
 * \see       seqan3::tuple_like
 */
//!\cond
template <typename tuple_t>
SEQAN3_CONCEPT tuple_size = requires (tuple_t v)
{
    { std::tuple_size<tuple_t>::value } -> size_t;
};
//!\endcond

/*!\interface seqan3::detail::tuple_get <>
 * \ingroup   core
 * \brief     Subconcept definition for seqan3::tuple_like to test for std::get-interface.
 * \see       seqan3::tuple_like
 */
//!\cond
template <typename tuple_t>
SEQAN3_CONCEPT tuple_get = requires (tuple_t & v, tuple_t const & v_c)
{
    requires std::tuple_size_v<tuple_t> > 0;

    typename std::tuple_element<0, tuple_t>::type;
    { get<0>(v)              } -> typename std::tuple_element<0, tuple_t>::type;
//     requires weakly_assignable_from<decltype(get<0>(v)), typename std::tuple_element<0, tuple_t>::type>;
    //TODO check that the previous returns something that can be assigned to
    // unfortunately std::assignable_from requires lvalue-reference, but we want to accept xvalues too (returned proxies)
    { get<0>(v_c)            } -> typename std::tuple_element<0, tuple_t>::type;
    { get<0>(std::move(v))   } -> typename std::tuple_element<0, tuple_t>::type;
    // TODO: The return type for std::tuple is wrong until gcc-8.0, for gcc > 8.0 this is fixed.
    { get<0>(std::move(v_c)) };// -> typename std::tuple_element<0, tuple_t>::type const &&;
};
//!\endcond

/*!\brief   Helper type trait function to check for std::totally_ordered on all elements of
 *          the given tuple type.
 * \ingroup core
 * \tparam  state_t   The last state of the fold operation.
 * \tparam  element_t The current processed element by the meta::fold operation.
 *
 * \returns std::true_type if std::totally_ordered<element_t> and state_t::value evaluate to `true`,
 *          std::false_type otherwise.
 */
template <typename state_t, typename element_t>
struct models_strict_totally_ordered
{
    //!\brief The resulting type definition.
    using type =  std::conditional_t<state_t::value && std::totally_ordered<element_t>,
                                    std::true_type,
                                    std::false_type>;
};

/*!\brief   Transformation trait to expose the tuple element types as seqan3::type_list
 * \ingroup core
 * \tparam  tuple_t The tuple to extract the element types from.
 *
 * \returns A seqan3::type_list over the element types of the given tuple.
 * \see seqan3::detail::tuple_type_list_t
 */
template <detail::tuple_size tuple_t>
struct tuple_type_list
{
protected:

    //!\brief Helper function to extract the types using the tuple elements.
    template <size_t ... Is>
    static constexpr auto invoke_to_type_list(std::index_sequence<Is...>)
    {
        return type_list<std::tuple_element_t<Is, tuple_t>...>{};
    }

public:
    //!\brief The generated seqan3::type_list.
    using type = decltype(invoke_to_type_list(std::make_index_sequence<std::tuple_size<tuple_t>::value>{}));
};

/*!\brief   Helper type for seqan3::detail::tuple_type_list
 * \ingroup core
 *
 * \see seqan3::detail::tuple_type_list
 */
template <detail::tuple_size tuple_t>
using tuple_type_list_t = typename tuple_type_list<tuple_t>::type;
} // namespace::seqan3

namespace seqan3
{

// ----------------------------------------------------------------------------
// tuple_like
// ----------------------------------------------------------------------------

/*!\interface   seqan3::tuple_like
 * \extends     std::totally_ordered
 * \ingroup     core
 * \brief       Whether a type behaves like a tuple.
 *
 * \details
 *
 * Types that meet this concept are for example std::tuple, std::pair, std::array, seqan3::pod_tuple, seqan3::record.
 * The std::totally_ordered will only be required if all types contained in the tuple like
 * data structure are them selfs strict totally ordered.
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
SEQAN3_CONCEPT tuple_like = detail::tuple_size<std::remove_reference_t<t>> && requires(t v)
{
    typename detail::tuple_type_list<remove_cvref_t<t>>::type;

    // NOTE(rrahn): To check the full tuple_concept including the get interface and the std::totally_ordered
    //              we need to make some assumptions. In general these checks can only be executed if the tuple is not
    //              empty. Furthermore, the std::totally_ordered can only be checked if all elements in the
    //              tuple are strict_totally_ordered. This is done, by the fold expression in the second part.
    requires (std::tuple_size<std::remove_reference_t<t>>::value == 0) ||
                detail::tuple_get<remove_cvref_t<t>> &&
                (!meta::fold<detail::tuple_type_list_t<remove_cvref_t<t>>,
                             std::true_type,
                             meta::quote_trait<detail::models_strict_totally_ordered>>::value ||
                std::totally_ordered<remove_cvref_t<t>>);
};
//!\endcond

} // namespace seqan3
