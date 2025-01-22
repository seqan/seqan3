// SPDX-FileCopyrightText: 2006-2025 Knut Reinert & Freie Universität Berlin
// SPDX-FileCopyrightText: 2016-2025 Knut Reinert & MPI für molekulare Genetik
// SPDX-License-Identifier: BSD-3-Clause

/*!\file
 * \author Hannes Hauswedell <hannes.hauswedell AT fu-berlin.de>
 * \author Lydia Buntrock <lydia.buntrock AT fu-berlin.de>
 * \brief Provides type traits for working with templates.
 */

#pragma once

#include <concepts>

#include <seqan3/utility/type_traits/detail/transformation_trait_or.hpp>

namespace seqan3::detail
{

// ----------------------------------------------------------------------------
// transfer_template_args_onto
// ----------------------------------------------------------------------------

//!\cond
template <typename source_type, template <typename...> typename target_template>
struct transfer_template_args_onto
{};
//!\endcond

/*!\brief Extracts a type template's **type** arguments and specialises another template with them.
 * \implements seqan3::transformation_trait
 * \ingroup core
 * \tparam source_template   The source type; must be a specialisation of a template.
 * \tparam target_template   The type template you wish to specialise.
 * \tparam source_arg_types  The **type** arguments to the source_template (deduced implicitly).
 * \see seqan3::detail::transfer_template_vargs_onto
 *
 * \details
 *
 * Among other use cases, it enables using the types contained in a seqan3::type_list to specialise another type
 * template.
 *
 * A type trait shortcut is also defined: seqan3::detail::transfer_template_args_onto_t
 *
 * This type trait works for templates that have **only type-arguments**. See
 * seqan3::detail::transfer_template_vargs_onto for a type trait that transfers non-type arguments. There is
 * no type trait that can handle a combination of type and non-type arguments.
 * If the `source_type` is a not a template class, e.g. an `int`, the member type `type` is not defined.
 *
 * ### Example
 *
 * \include test/snippet/core/detail/template_inspection_usage.cpp
 */
template <template <typename...> typename source_template,
          template <typename...>
          typename target_template,
          typename... source_arg_types>
    requires requires () { typename target_template<source_arg_types...>; }
struct transfer_template_args_onto<source_template<source_arg_types...>, target_template>
{
    //!\brief The return type: the target type specialised by the unpacked types in the list.
    using type = target_template<source_arg_types...>;
};

/*!\brief Shortcut for seqan3::detail::transfer_template_args_onto (transformation_trait shortcut).
 * \ingroup core
 * \see seqan3::detail::transfer_template_args_onto
 */
template <typename source_type, template <typename...> typename target_template>
using transfer_template_args_onto_t = typename transfer_template_args_onto<source_type, target_template>::type;

// ----------------------------------------------------------------------------
// transfer_template_vargs_onto
// ----------------------------------------------------------------------------

//!\cond
template <typename source_type, template <auto...> typename target_template>
struct transfer_template_vargs_onto
{};
//!\endcond

/*!\brief Extracts a type template's **non-type** arguments and specialises another template with them.
 * \implements seqan3::transformation_trait
 * \ingroup core
 * \tparam source_template   The source type; must be a specialisation of a template.
 * \tparam target_template   The type template you wish to specialise.
 * \tparam source_varg_types The **non-type** arguments to the source_template (deduced implicitly).
 * \see seqan3::detail::transfer_template_vargs_onto
 *
 * \details
 *
 * A shortcut is also defined: seqan3::detail::transfer_template_vargs_onto_t
 *
 * This transformation trait works for templates that have **only non-type-arguments**. See
 * seqan3::detail::transfer_template_args_onto for a transformation trait that transfers type arguments. There is
 * no transformation trait that can handle a combination of type and non-type arguments.
 * If the `source_type` is a not a template class, e.g. an `int`, the member type `type` is not defined.
 */
template <template <auto...> typename source_template,
          template <auto...>
          typename target_template,
          auto... source_varg_types>
    requires requires () { typename target_template<source_varg_types...>; }
struct transfer_template_vargs_onto<source_template<source_varg_types...>, target_template>
{
    //!\brief The return type: the target type specialised by the unpacked types in the list.
    using type = target_template<source_varg_types...>;
};

/*!\brief Shortcut for seqan3::detail::transfer_template_vargs_onto (transformation_trait shortcut).
 * \ingroup core
 * \see seqan3::detail::transfer_template_vargs_onto
 */
template <typename source_type, template <auto...> typename target_template>
using transfer_template_vargs_onto_t = typename transfer_template_vargs_onto<source_type, target_template>::type;

// ----------------------------------------------------------------------------
// is_type_specialisation_of_v
// ----------------------------------------------------------------------------

/*!\brief Determines whether a source_type is a specialisation of another template.
 * \implements seqan3::unary_type_trait
 * \ingroup core
 * \tparam source_type      The source type.
 * \tparam target_template  The type template you wish to compare against (must take only types as template arguments).
 *
 * \details
 *
 * ### Example
 *
 * \include test/snippet/core/detail/template_inspection_usage_2.cpp
 */
template <typename source_t, template <typename...> typename target_template>
struct is_type_specialisation_of : public std::false_type
{};

//!\overload
template <template <typename...> typename source_t, typename... source_args>
struct is_type_specialisation_of<source_t<source_args...>, source_t> : public std::true_type
{};

/*!\brief Helper variable template for seqan3::detail::is_type_specialisation_of (unary_type_trait shortcut).
 * \relates seqan3::detail::is_type_specialisation_of
 * \tparam source_type      The source type.
 * \tparam target_template  The type template you wish to compare against (must take only types as template arguments).
 */
template <typename source_t, template <typename...> typename target_template>
inline constexpr bool is_type_specialisation_of_v = is_type_specialisation_of<source_t, target_template>::value;

// ----------------------------------------------------------------------------
// is_value_specialisation_of_v
// ----------------------------------------------------------------------------

//!\cond
template <typename source_t, template <auto...> typename target_template>
struct is_value_specialisation_of : std::false_type
{};
//!\endcond

/*!\brief Determines whether a source_type is a specialisation of another template.
 * \implements seqan3::unary_type_trait
 * \ingroup core
 * \tparam source_type      The source type.
 * \tparam target_template  The type template you wish to compare against (must take only non-types as template
 * arguments).
 * \see seqan3::detail::is_type_specialisation_of
 * \see seqan3::detail::is_value_specialisation_of_v
 */
template <template <auto...> typename source_t, auto... source_args>
struct is_value_specialisation_of<source_t<source_args...>, source_t> : public std::true_type
{};

/*!\brief Helper variable template for seqan3::detail::is_value_specialisation_of (unary_type_trait shortcut).
 * \relates seqan3::detail::is_value_specialisation_of
 * \tparam source_type      The source type.
 * \tparam target_template  The type template you wish to compare against (must take only types as template arguments).
 */
template <typename source_t, template <auto...> typename target_template>
inline constexpr bool is_value_specialisation_of_v = is_value_specialisation_of<source_t, target_template>::value;

/*!\brief Exposes `templ_t<spec_t...>` if that is valid, otherwise `fallback_t`.
 * \implements seqan3::transformation_trait
 * \see seqan3::detail::valid_template_spec_or_t
 * \ingroup core
 * \tparam fallback_t The fallback type.
 * \tparam templ_t    The type template that should be specialised.
 * \tparam spec_t     The specialisation for the type template.
 */
template <typename fallback_t, template <typename...> typename templ_t, typename... spec_t>
struct valid_template_spec_or
{
    //!\brief The resulting type.
    using type = fallback_t;
};

//!\overload
template <typename fallback_t, template <typename...> typename templ_t, typename... spec_t>
    requires requires { typename templ_t<spec_t...>; }
struct valid_template_spec_or<fallback_t, templ_t, spec_t...>
{
    //!\brief The resulting type.
    using type = templ_t<spec_t...>;
};

/*!\brief Helper for seqan3::detail::valid_template_spec_or (transformation_trait shortcut).
 * \ingroup core
 * \see seqan3::detail::valid_template_spec_or
 * \tparam fallback_t The fallback type.
 * \tparam templ_t    The type template that should be specialised.
 * \tparam spec_t     The specialisation for the type template.
 */
template <typename fallback_t, template <typename...> typename templ_t, typename... spec_t>
using valid_template_spec_or_t = typename valid_template_spec_or<fallback_t, templ_t, spec_t...>::type;

// ----------------------------------------------------------------------------
// template_specialisation_of
// ----------------------------------------------------------------------------

/*!\interface seqan3::detail::template_specialisation_of <>
 * \brief Provides concept `seqan3::template_specialisation_of<mytype, [...]>` for checking the type specialisation of
 *        some type with a given template, for example a specialized `type_list<float>` with the `type_list` template.
 * \ingroup core
 *
 * \tparam mytype           The query type.
 * \tparam type_template    The type template you wish to compare against mytype.
 *
 * \see seqan3::detail::is_type_specialisation_of_v
 *
 * ### Example
 *
 * \include test/snippet/core/detail/template_inspection_usage_3.cpp
 */
//!\cond
template <typename mytype, template <typename...> typename type_template>
concept template_specialisation_of = is_type_specialisation_of_v<mytype, type_template>;

//!\endcond

// ----------------------------------------------------------------------------
// strip_type_identity
// ----------------------------------------------------------------------------

/*!\brief A transformation trait shortcut that returns the type inside a std::type_identity or the type itself.
 * \ingroup core
 * \tparam t The type to operate on.
 */
template <typename t>
using strip_type_identity_t =
    std::conditional_t<is_type_specialisation_of_v<t, std::type_identity>, transformation_trait_or_t<t, void>, t>;
} // namespace seqan3::detail
