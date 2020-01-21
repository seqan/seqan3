// -----------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2020, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2020, Knut Reinert & MPI für molekulare Genetik
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at: https://github.com/seqan/seqan3/blob/master/LICENSE.md
// -----------------------------------------------------------------------------------------------------

/*!\file
 * \brief Provides seqan3::search_cfg::parallel configuration.
 * \author Svenja Mehringer <svenja.mehringer AT fu-berlin.de>
 */

#pragma once

#include <seqan3/core/algorithm/configuration_element_parallel_mode.hpp>
#include <seqan3/search/configuration/detail.hpp>

namespace seqan3::search_cfg
{
/*!\brief Enables the parallel execution of the search algorithm if possible for the given configuration.
 * \ingroup search_configuration
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
using parallel = seqan3::detail::parallel_mode<std::integral_constant<detail::search_config_id,
                                                                      detail::search_config_id::parallel>>;

} // namespace seqan3::search_cfg
