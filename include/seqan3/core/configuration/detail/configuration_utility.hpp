// -----------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2020, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2020, Knut Reinert & MPI für molekulare Genetik
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at: https://github.com/seqan/seqan3/blob/master/LICENSE.md
// -----------------------------------------------------------------------------------------------------

/*!\file
 * \brief Provides various auxiliary functions with which parts of the configurations can be checked.
 * \author Rene Rahn <rene.rahn AT fu-berlin.de>
 */

#pragma once

#include <seqan3/core/detail/template_inspection.hpp>

namespace seqan3::detail
{

// ----------------------------------------------------------------------------
// Metafunction is_same_configuration_f
// ----------------------------------------------------------------------------

/*!\brief Helper meta function to check if a template type is contained in a seqan3::configuration.
 * \ingroup algorithm
 *
 * \details
 *
 * This helper meta function is used to provide the `get` and `get_or` interface for template template types.
 */
template <template <typename ...> typename query_t>
struct is_same_configuration_f
{
    /*!\brief A type template that evaluates to std::true_type if the given type is a specialization of `query_t`,
     *        otherwise std::false_type.
     * \tparam compare_type The type to compare against `query_t`.
     */
    template <typename compare_type>
    using invoke = is_type_specialisation_of<compare_type, query_t>;
};

} // namespace seqan3::detail
