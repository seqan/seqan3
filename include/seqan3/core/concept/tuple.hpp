// ============================================================================
//                 SeqAn - The Library for Sequence Analysis
// ============================================================================
//
// Copyright (c) 2006-2018, Knut Reinert & Freie Universitaet Berlin
// Copyright (c) 2016-2018, Knut Reinert & MPI Molekulare Genetik
// All rights reserved.
//
// Redistribution and use in source and binary forms, with or without
// modification, are permitted provided that the following conditions are met:
//
//     * Redistributions of source code must retain the above copyright
//       notice, this list of conditions and the following disclaimer.
//     * Redistributions in binary form must reproduce the above copyright
//       notice, this list of conditions and the following disclaimer in the
//       documentation and/or other materials provided with the distribution.
//     * Neither the name of Knut Reinert or the FU Berlin nor the names of
//       its contributors may be used to endorse or promote products derived
//       from this software without specific prior written permission.
//
// THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
// AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
// IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
// ARE DISCLAIMED. IN NO EVENT SHALL KNUT REINERT OR THE FU BERLIN BE LIABLE
// FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL
// DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR
// SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER
// CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT
// LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY
// OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH
// DAMAGE.
//
// ============================================================================

/*!\file
 * \brief Provides seqan3::tuple_like_concept.
 * \author Hannes Hauswedell <hannes.hauswedell AT fu-berlin.de>
 */

#pragma once

#include <tuple>
#include <type_traits>

#include <seqan3/core/metafunction/basic.hpp>
#include <seqan3/core/metafunction/template_inspection.hpp>
#include <seqan3/core/type_list.hpp>
#include <seqan3/core/pod_tuple.hpp>
#include <seqan3/std/concepts>

namespace seqan3::detail
{

/*!\interface seqan3::detail::tuple_size_concept <>
 * \ingroup   core
 * \brief     Subconcept definition for seqan3::tuple_like_concept to test for std::tuple_size-interface.
 * \see       seqan3::tuple_like_concept
 */
//!\cond
template <typename tuple_t>
concept tuple_size_concept = requires (tuple_t v)
{
    { std::tuple_size<tuple_t>::value } -> size_t;
};
//!\endcond

/*!\interface seqan3::detail::tuple_get_concept <>
 * \ingroup   core
 * \brief     Subconcept definition for seqan3::tuple_like_concept to test for std::get-interface.
 * \see       seqan3::tuple_like_concept
 */
//!\cond
template <typename tuple_t>
concept tuple_get_concept = requires (tuple_t & v, tuple_t const & v_c)
{
    requires std::tuple_size_v<tuple_t> > 0;

    typename std::tuple_element<0, tuple_t>::type;
    { std::get<0>(v)              } -> typename std::tuple_element<0, tuple_t>::type &;
    { std::get<0>(v_c)            } -> typename std::tuple_element<0, tuple_t>::type const &;
    { std::get<0>(std::move(v))   } -> typename std::tuple_element<0, tuple_t>::type &&;
    // TODO: The return type for std::tuple is wrong until gcc-8.0, for gcc > 8.0 this is fixed.
    { std::get<0>(std::move(v_c)) };// -> typename std::tuple_element<0, tuple_t>::type const &&;
};
//!\endcond

/*!\brief   Helper type trait function to check for std::StrictTotallyOrdered on all elements of
 *          the given tuple type.
 * \ingroup core
 * \tparam  state_t   The last state of the fold operation.
 * \tparam  element_t The current processed element by the meta::fold operation.
 *
 * \returns std::true_type if std::StrictTotallyOrdered<element_t> and state_t::value evaluate to `true`,
 *          std::false_type otherwise.
 */
template <typename state_t, typename element_t>
struct models_strict_totally_ordered
{
    //!\brief The resulting type definition.
    using type =  std::conditional_t<state_t::value && std::StrictTotallyOrdered<element_t>,
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
template <detail::tuple_size_concept tuple_t>
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
template <detail::tuple_size_concept tuple_t>
using tuple_type_list_t = typename tuple_type_list<tuple_t>::type;
} // namespace::seqan3

namespace seqan3
{

// ----------------------------------------------------------------------------
// tuple_like_concept
// ----------------------------------------------------------------------------

/*!\interface   seqan3::tuple_like_concept
 * \extends     std::StrictTotallyOrdered
 * \ingroup     core
 * \brief       Whether a type behaves like a tuple.
 *
 * \details
 *
 * Types that meet this concept are for example std::tuple, std::pair, std::array, seqan3::pod_tuple, seqan3::record.
 * The std::StrictTotallyOrdered will only be required if all types contained in the tuple like
 * data structure are them selfs strict totally ordered.
 */
/*!\name Requirements for seqan3::tuple_like_concept
 * \brief You can expect these (meta-)functions on all types that implement seqan3::tuple_like_concept.
 * \{
 */
/*!\var         size_t std::tuple_size_v<type>
 * \brief       A unary type trait that holds the number of elements in the tuple.
 * \tparam      type The tuple-like type.
 * \relates     seqan3::tuple_like_concept
 *
 * \details
 * \attention This is a concept requirement, not an actual function (however types satisfying this concept
 * will provide an implementation).
 */
/*!\typedef     std::tuple_elment_t<i, type>
 * \brief       A transformation trait that holds the type of elements in the tuple.
 * \tparam      i Index of the queried element type.
 * \tparam      type The tuple-like type.
 * \relates     seqan3::tuple_like_concept
 *
 * \details
 * \attention This is a concept requirement, not an actual function (however types satisfying this concept
 * will provide an implementation).
 * \attention This constraint is **not enforced** since empty tuples are valid.
 */
/*!\fn              auto && std::get<i>(type && val)
 * \brief           Return the i-th element of the tuple.
 * \relates         seqan3::tuple_like_concept
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
concept tuple_like_concept = detail::tuple_size_concept<std::remove_reference_t<t>> && requires(t v)
{
    typename detail::tuple_type_list<remove_cvref_t<t>>::type;

    // NOTE(rrahn): To check the full tuple_concept including the get interface and the std::StrictTotallyOrdered
    //              we need to make some assumptions. In general these checks can only be executed if the tuple is not
    //              empty. Furthermore, the std::StrictTotallyOrdered can only be checked if all elements in the
    //              tuple are strict_totally_ordered. This is done, by the fold expression in the second part.
    requires (std::tuple_size<std::remove_reference_t<t>>::value == 0) ||
                detail::tuple_get_concept<remove_cvref_t<t>> &&
                (!meta::fold<detail::tuple_type_list_t<remove_cvref_t<t>>,
                             std::true_type,
                             meta::quote_trait<detail::models_strict_totally_ordered>>::value ||
                std::StrictTotallyOrdered<remove_cvref_t<t>>);
};
//!\endcond

} // namespace seqan3
