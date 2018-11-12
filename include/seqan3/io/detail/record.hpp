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

/*!\cond DEV
 * \file
 * \brief Provides auxiliary data structures and functions for seqan3::record and seqan3::fields.
 * \author Hannes Hauswedell <hannes.hauswedell AT fu-berlin.de>
 * \endcond
 */

#pragma once

#include <range/v3/view/repeat.hpp>

#include <seqan3/core/metafunction/basic.hpp>
#include <seqan3/io/record.hpp>
#include <seqan3/std/ranges>

namespace seqan3::detail
{

// ----------------------------------------------------------------------------
// Fields 
// ----------------------------------------------------------------------------

/*!\brief Auxiliary concept that checks whether a type is a specialisation of seqan3::fields.
 * \ingroup io
 * \relates seqan3::fields
 */
template <typename t>
concept Fields  = is_value_specialisation_of_v<t, fields>;

// ----------------------------------------------------------------------------
// select_types_with_ids
// ----------------------------------------------------------------------------

/*!\brief A type metafunction that selects a subset of types based on their IDs.
 * \ingroup io
 * \tparam field_types          The types of the fields available to the record in a seqan3::type_list.
 * \tparam field_types_as_ids   A seqan3::fields type with seqan3::field IDs corresponding to field_types.
 * \tparam selected_field_ids   A seqan3::fields type with the subset (and order) of the fields selected.
 * \tparam field_no             The field we are currently processing (defaults to 0).
 * \tparam return_types         The type pack being aggregated (empty at start).
 *
 * Given a list of types and corresponding IDs; and given a selection (and possibly different order) of IDs, return
 * the types corresponding to that selection and in that order.
 *
 * This metafunction recurses over the `selected_field_ids` and retrieves the corresponding typenames from
 * `field_types` via their identifer in `field_types_as_ids`. It recursively builds up `return_types` which it packs
 * into a seqan3::type_list once the end of selected_field_ids is reached.
 *
 * ### Example
 *
 * ```cpp
 * using types         = type_list<std::string, dna4_vector, std::vector<phred42>>;
 * using types_as_ids  = fields<field::ID,      field::SEQ,  field::QUAL>;
 * using selected_ids  = fields<field::QUAL, field::ID>;
 *
 * using selected_types = detail::select_types_with_ids_t<types, types_as_ids, selected_ids>;
 * // resolves to type_list<std::vector<phred42>, std::string>>
 * ```
 */
template <typename field_types,
          typename field_types_as_ids,
          typename selected_field_ids,
          size_t field_no = 0,
          typename ... return_types>
struct select_types_with_ids                               // unconstrained template is recursion anchor
{
    //!\brief The return type.
    using type = type_list<return_types...>;
};

/*!\brief Type metafunction shortcut for seqan3::select_types_with_ids.
 * \ingroup io
 * \relates select_types_with_ids
 */
template <typename field_types,
          typename field_types_as_ids,
          typename selected_field_ids,
          size_t field_no = 0,
          typename ... return_types>
using select_types_with_ids_t = typename select_types_with_ids<field_types,
                                                               field_types_as_ids,
                                                               selected_field_ids,
                                                               field_no,
                                                               return_types...>::type;
//!\cond
template <typename field_types,
          typename field_types_as_ids,
          typename selected_field_ids,
          size_t field_no,
          typename ... return_types>
    requires field_no < selected_field_ids::as_array.size() // perform recursion while not at end
struct select_types_with_ids<field_types, field_types_as_ids, selected_field_ids, field_no, return_types...>
{
    static_assert(field_types_as_ids::contains(selected_field_ids::as_array[field_no]),
                  "You selected a field that was not in field_types_as_ids.");

    // call this metafunction again, but increase index by one and append a type to the returned type list.
    using type = select_types_with_ids_t<field_types,
                                         field_types_as_ids,
                                         selected_field_ids,
                                         field_no + 1,
                                         return_types ...,
                                         meta::at_c<field_types,
                                                    field_types_as_ids::index_of(selected_field_ids::as_array[field_no])>>;

};
//!\endcond


// ----------------------------------------------------------------------------
// get_or_ignore
// ----------------------------------------------------------------------------

/*!\addtogroup io
 *!\{
 */
//!\brief Access an element in a std::tuple or seqan3::record; return reference to std::ignore if not contained.
template <field f, typename field_types, typename field_ids>
auto & get_or_ignore(record<field_types, field_ids> & r)
{
    if constexpr (field_ids::contains(f))
        return std::get<field_ids::index_of(f)>(r);
    else
        return std::ignore;
}

//!\copydoc seqan3::detail::get_or_ignore
template <field f, typename field_types, typename field_ids>
auto const & get_or_ignore(record<field_types, field_ids> const & r)
{
    if constexpr (field_ids::contains(f))
        return std::get<field_ids::index_of(f)>(r);
    else
        return std::ignore;
}

//!\copydoc seqan3::detail::get_or_ignore
template <size_t i, typename ... types>
auto & get_or_ignore(std::tuple<types...> & t)
{
    if constexpr (i < sizeof...(types))
        return std::get<i>(t);
    else
        return std::ignore;
}

//!\copydoc seqan3::detail::get_or_ignore
template <size_t i, typename ... types>
auto const & get_or_ignore(std::tuple<types...> const & t)
{
    if constexpr (i < sizeof...(types))
        return std::get<i>(t);
    else
        return std::ignore;
}
//!\}


// ----------------------------------------------------------------------------
// get_or
// ----------------------------------------------------------------------------

/*!\addtogroup io
 *!\{
 */
//!\brief Access an element in a std::tuple or seqan3::record; return or_value if not contained.
template <field f, typename field_types, typename field_ids, typename or_type>
decltype(auto) get_or(record<field_types, field_ids> & r, or_type && or_value)
{
    if constexpr (field_ids::contains(f))
        return std::get<field_ids::index_of(f)>(r);
    else
        return std::forward<or_type>(or_value);
}

//!\copydoc seqan3::detail::get_or
template <field f, typename field_types, typename field_ids, typename or_type>
decltype(auto) get_or(record<field_types, field_ids> const & r, or_type && or_value)
{
    if constexpr (field_ids::contains(f))
        return std::get<field_ids::index_of(f)>(r);
    else
        return std::forward<or_type>(or_value);
}

//!\copydoc seqan3::detail::get_or
template <size_t i, typename or_type, typename ... types>
decltype(auto) get_or(std::tuple<types...> & t, or_type && or_value)
{
    if constexpr (i < sizeof...(types))
        return std::get<i>(t);
    else
        return std::forward<or_type>(or_value);
}

//!\copydoc seqan3::detail::get_or
template <size_t i, typename or_type, typename ... types>
decltype(auto) get_or(std::tuple<types...> const & t, or_type && or_value)
{
    if constexpr (i < sizeof...(types))
        return std::get<i>(t);
    else
        return std::forward<or_type>(or_value);
}
//!\}

// ----------------------------------------------------------------------------
// range_wrap_ignore
// ----------------------------------------------------------------------------

//!\brief Pass through the reference to the argument in case the argument satisfies std::ranges::InputRange.
template <std::ranges::InputRange rng_t>
inline auto & range_wrap_ignore(rng_t & range)
{
    return range;
}

/*!\brief If the argument is std::ignore, return an infinite range of std::ignore values.
 * \details
 *
 * This function can be used in combination with seqan3::detail::get_or_ignore to ensure same dimensionality of
 * the returned type, even for fields not present in the record / tuple.
 */
inline auto range_wrap_ignore(ignore_t const &)
{
    return ranges::view::repeat(std::ignore);
}

} // namespace seqan3::detail
