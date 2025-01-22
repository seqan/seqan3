// SPDX-FileCopyrightText: 2006-2025 Knut Reinert & Freie Universität Berlin
// SPDX-FileCopyrightText: 2016-2025 Knut Reinert & MPI für molekulare Genetik
// SPDX-License-Identifier: BSD-3-Clause

/*!\file
 * \brief Provides seqan3::search_cfg::parallel configuration.
 * \author Svenja Mehringer <svenja.mehringer AT fu-berlin.de>
 */

#pragma once

#include <seqan3/core/configuration/detail/configuration_element_parallel_mode.hpp>
#include <seqan3/search/configuration/detail.hpp>

namespace seqan3::search_cfg
{
/*!\brief Enables the parallel execution of the search algorithm if possible for the given configuration.
 * \ingroup search_configuration
 * \see search_configuration
 *
 * \details
 *
 * With this configuration you can enable the parallel execution of the search algorithm.
 *
 * The config element takes the number of threads as a parameter, which must be greater than `0`.
 *
 * ### Example
 *
 * \include test/snippet/search/configuration_parallel.cpp
 */
using parallel = seqan3::detail::parallel_mode<
    std::integral_constant<seqan3::detail::search_config_id, seqan3::detail::search_config_id::parallel>>;

} // namespace seqan3::search_cfg
