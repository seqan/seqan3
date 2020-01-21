// -----------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2020, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2020, Knut Reinert & MPI für molekulare Genetik
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at: https://github.com/seqan/seqan3/blob/master/LICENSE.md
// -----------------------------------------------------------------------------------------------------

/*!\file
 * \brief Provides seqan3::align_cfg::parallel configuration.
 * \author Lydia Buntrock <lydia.buntrock AT fu-berlin.de>
 * \author Jörg Winkler <j.winkler AT fu-berlin.de>
 * \author Svenja Mehringer <svenja.mehringer AT fu-berlin.de>
 */

#pragma once

#include <seqan3/alignment/configuration/detail.hpp>
#include <seqan3/core/algorithm/configuration_element_parallel_mode.hpp>

namespace seqan3::align_cfg
{
/*!\brief Enables the parallel execution of the alignment algorithm if possible for the given configuration.
 * \ingroup alignment_configuration
 *
 * \details
 *
 * With this configuration you can enable the parallel execution of the pairwise sequence alignment. This means that for
 * a batch of pairwise sequence alignments the specified number of threads will be spawned to compute them in parallel.
 * Note that only independent alignment computations can be executed in parallel, i.e. you use this method when
 * computing a batch of alignments rather than executing them separately.
 * Depending on your processor architecture you can gain a significant speed-up.
 *
 * The value represents the number of threads to be used and must be greater than `0`.
 *
 * ### Example
 *
 * \include test/snippet/alignment/configuration/align_cfg_parallel_example.cpp
 */
using parallel = seqan3::detail::parallel_mode<std::integral_constant<detail::align_config_id,
                                                                      detail::align_config_id::parallel>>;

} // namespace seqan3::align_cfg
