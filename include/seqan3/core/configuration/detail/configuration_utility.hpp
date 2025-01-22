// SPDX-FileCopyrightText: 2006-2025 Knut Reinert & Freie Universität Berlin
// SPDX-FileCopyrightText: 2016-2025 Knut Reinert & MPI für molekulare Genetik
// SPDX-License-Identifier: BSD-3-Clause

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
 * \ingroup core_configuration
 *
 * \details
 *
 * This helper meta function is used to provide the `get` and `get_or` interface for template template types.
 */
template <template <typename...> typename query_t>
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
