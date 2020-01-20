// -----------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2020, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2020, Knut Reinert & MPI für molekulare Genetik
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at: https://github.com/seqan/seqan3/blob/master/LICENSE.md
// -----------------------------------------------------------------------------------------------------

/*!\cond DEV
 * \file
 * \brief Provides auxiliary data structures and functions for seqan3::record and seqan3::fields.
 * \author Hannes Hauswedell <hannes.hauswedell AT fu-berlin.de>
 * \endcond
 */

#pragma once

#include <seqan3/core/concept/tuple.hpp>
#include <seqan3/core/type_list/traits.hpp>
#include <seqan3/io/record.hpp>
#include <seqan3/range/views/repeat.hpp>
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
SEQAN3_CONCEPT fields_specialisation = is_value_specialisation_of_v<t, fields>;

// ----------------------------------------------------------------------------
// select_types_with_ids
// ----------------------------------------------------------------------------

/*!\brief Exposes a subset of types as a seqan3::type_list selected based on their IDs.
 * \implements seqan3::transformation_trait
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
 * This transformation trait recurses over the `selected_field_ids` and retrieves the corresponding typenames from
 * `field_types` via their identifer in `field_types_as_ids`. It recursively builds up `return_types` which it packs
 * into a seqan3::type_list once the end of selected_field_ids is reached.
 *
 * ### Example
 *
 * \include test/snippet/io/detail/detail_record.cpp
 */
template <typename field_types,
          typename field_types_as_ids,
          typename selected_field_ids,
          size_t field_no = 0,
          typename ...return_types>
struct select_types_with_ids                               // unconstrained template is recursion anchor
{
    //!\brief The return type.
    using type = type_list<return_types...>;
};

/*!\brief Shortcut for seqan3::select_types_with_ids (transformation_trait shortcut).
 * \ingroup io
 * \relates select_types_with_ids
 */
template <typename field_types,
          typename field_types_as_ids,
          typename selected_field_ids,
          size_t field_no = 0,
          typename ...return_types>
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
          typename ...return_types>
    requires field_no < selected_field_ids::as_array.size() // perform recursion while not at end
struct select_types_with_ids<field_types, field_types_as_ids, selected_field_ids, field_no, return_types...>
{
    static_assert(field_types_as_ids::contains(selected_field_ids::as_array[field_no]),
                  "You selected a field that was not in field_types_as_ids.");

    // call this type trait again, but increase index by one and append a type to the returned type list.
    using type = select_types_with_ids_t<field_types,
                                         field_types_as_ids,
                                         selected_field_ids,
                                         field_no + 1,
                                         return_types ...,
                                         list_traits::at<field_types_as_ids::index_of(
                                             selected_field_ids::as_array[field_no]), field_types>>;

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
template <size_t i, template <tuple_like ...types_> typename tuple_like_t, typename ...types>
auto & get_or_ignore(tuple_like_t<types...> & t)
{
    if constexpr (i < sizeof...(types))
        return std::get<i>(t);
    else
        return std::ignore;
}

//!\copydoc seqan3::detail::get_or_ignore
template <size_t i, template <tuple_like ...types_> typename tuple_like_t, typename ...types>
auto const & get_or_ignore(tuple_like_t<types...> const & t)
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
template <size_t i, typename or_type, typename ...types>
decltype(auto) get_or(std::tuple<types...> & t, or_type && or_value)
{
    if constexpr (i < sizeof...(types))
        return std::get<i>(t);
    else
        return std::forward<or_type>(or_value);
}

//!\copydoc seqan3::detail::get_or
template <size_t i, typename or_type, typename ...types>
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

//!\brief Pass through the reference to the argument in case the argument satisfies std::ranges::input_range.
template <std::ranges::input_range rng_t>
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
    return views::repeat(std::ignore);
}

} // namespace seqan3::detail
