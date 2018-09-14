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
 * \author Hannes Hauswedell <hannes.hauswedell AT fu-berlin.de>
 * \brief Provides seqan3::type_list and auxiliary metafunctions.
 */

#pragma once

#include <meta/meta.hpp>

#include <seqan3/core/platform.hpp>

namespace seqan3::detail
{

// ----------------------------------------------------------------------------
// transfer_template_args_onto
// ----------------------------------------------------------------------------

/*!\brief Type metafunction that extracts a type template's **types** arguments and specialises another template
 * with them [metafunction declaration].
 * \tparam source_type      The source type.
 * \tparam target_template  The type template you wish to specialise.
 * \ingroup metafunction
 * \see seqan3::detail::transfer_template_vargs_onto
 *
 * Among other use cases, it enables using the types contained in a seqan3::type_list to specialise another type
 * template.
 *
 * A metafunction shortcut is also defined: seqan3::detail::transfer_template_args_onto_t
 *
 * This metafunction works for templates that have **only type-arguments**. See
 * seqan3::detail::transfer_template_vargs_onto for a metafunction that transfers non-type arguments. There is
 * no metafunction that can handle a combination of type and non-type arguments.
 * If the `source_type` is a not a template class, e.g. an `int`, the return type defaults to `void`.
 *
 * ### Example
 *
 * \snippet test/snippet/core/metafunction/template_inspection.cpp usage
 */
template <typename source_type, template <typename ...> typename target_template>
struct transfer_template_args_onto
{
    //!\brief The return type: set to void for non template types.
    using type = void;
};

/*!\brief Type metafunction that extracts a type template's **type** arguments and specialises another template
 * with them [metafunction definition].
 * \tparam source_template   The source type; must be a specialisation of a template.
 * \tparam target_template   The type template you wish to specialise.
 * \tparam source_arg_types  The **type** arguments to the source_template (deduced implicitly).
 * \ingroup metafunction
 * \see seqan3::detail::transfer_template_args_onto
 *
 * The only viable definition assumes a source type that is a fully specialised template with only **type**
 * arguments.
 */
template <template <typename ...> typename source_template,
          template <typename ...> typename target_template,
          typename ... source_arg_types>
struct transfer_template_args_onto<source_template<source_arg_types...>, target_template>
{
    //!\brief The return type: the target type specialised by the unpacked types in the list.
    using type = target_template<source_arg_types...>;
};

/*!\brief Type metafunction shortcut for seqan3::detail::transfer_template_args_onto.
 * \ingroup metafunction
 * \see seqan3::detail::transfer_template_args_onto
 */
template <typename source_type, template <typename ...> typename target_template>
using transfer_template_args_onto_t = typename transfer_template_args_onto<source_type, target_template>::type;

// ----------------------------------------------------------------------------
// transfer_template_vargs_onto
// ----------------------------------------------------------------------------

/*!\brief Type metafunction that extracts a type template's **non-type** arguments and specialises another template
 * with them [metafunction declaration].
 * \tparam source_type      The source type.
 * \tparam target_template  The type template you wish to specialise.
 * \ingroup metafunction
 * \see seqan3::detail::transfer_template_args_onto
 *
 * There is a shortcut for this metafunction: seqan3::detail::transfer_template_args_onto_t
 *
 * This metafunction works for templates that have **only non-type-arguments**. See
 * seqan3::detail::transfer_template_args_onto for a metafunction that transfers type arguments. There is
 * no metafunction that can handle a combination of type and non-type arguments.
 * If the `source_type` is a not a template class, e.g. an `int`, the return type defaults to `void`.
 */

template <typename source_type, template <auto ...> typename target_template>
struct transfer_template_vargs_onto
{
    //!\brief The return type: set to void for non template types.
    using type = void;
};

/*!\brief Type metafunction that extracts a type template's **non-type** arguments and specialises another template
 * with them [metafunction definition].
 * \tparam source_template   The source type; must be a specialisation of a template.
 * \tparam target_template   The type template you wish to specialise.
 * \tparam source_varg_types The **non-type** arguments to the source_template (deduced implicitly).
 * \ingroup metafunction
 * \see seqan3::detail::transfer_template_vargs_onto
 *
 * The only viable definition assumes a source type that is a fully specialised template with only **non-type**
 * arguments.
 */
template <template <auto ...> typename source_template,
          template <auto ...> typename target_template,
          auto ... source_varg_types>
struct transfer_template_vargs_onto<source_template<source_varg_types...>, target_template>
{
    //!\brief The return type: the target type specialised by the unpacked types in the list.
    using type = target_template<source_varg_types...>;
};

/*!\brief Type metafunction shortcut for seqan3::detail::transfer_template_vargs_onto.
 * \ingroup metafunction
 * \see seqan3::detail::transfer_template_vargs_onto
 */
template <typename source_type, template <auto ...> typename target_template>
using transfer_template_vargs_onto_t = typename transfer_template_vargs_onto<source_type, target_template>::type;

// ----------------------------------------------------------------------------
// is_type_specialisation_of_v
// ----------------------------------------------------------------------------

/*!\brief Value metafunction that returns whether a source_type is a specialisation of another template.
 * \tparam source_type      The source type.
 * \tparam target_template  The type template you wish to compare against (must take only types as template arguments).
 * \ingroup metafunction
 * \see seqan3::detail::is_value_specialisation_of
 * \see seqan3::detail::is_type_specialisation_of_v
 *
 * \details
 *
 * ### Example
 *
 * \snippet test/snippet/core/metafunction/template_inspection.cpp usage_2
 */
template <typename source_t, template <typename ...> typename target_template>
struct is_type_specialisation_of :
    std::is_same<source_t, transfer_template_args_onto_t<source_t, target_template>>
{};

/*!\brief Helper variable template for seqan3::detail::is_type_specialisation_of.
 * \tparam source_type      The source type.
 * \tparam target_template  The type template you wish to compare against (must take only types as template arguments).
 * \ingroup metafunction
 */
template <typename source_t, template <typename ...> typename target_template>
inline constexpr bool is_type_specialisation_of_v = is_type_specialisation_of<source_t, target_template>::value;

// ----------------------------------------------------------------------------
// is_value_specialisation_of_v
// ----------------------------------------------------------------------------

/*!\brief Value metafunction that returns whether a source_type is a specialisation of another template.
 * \tparam source_type      The source type.
 * \tparam target_template  The type template you wish to compare against (must take only non-types as template
 * arguments).
 * \ingroup metafunction
 * \see seqan3::detail::is_type_specialisation_of
 * \see seqan3::detail::is_value_specialisation_of_v
 */
template <typename source_t, template <auto ...> typename target_template>
struct is_value_specialisation_of :
    std::is_same<source_t, transfer_template_vargs_onto_t<source_t, target_template>>
{};

/*!\brief Helper variable template for seqan3::detail::is_value_specialisation_of.
 * \tparam source_type      The source type.
 * \tparam target_template  The type template you wish to compare against (must take only types as template arguments).
 * \ingroup metafunction
 */
template <typename source_t, template <auto ...> typename target_template>
inline constexpr bool is_value_specialisation_of_v = is_value_specialisation_of<source_t, target_template>::value;

} // namespace seqan3::detail
