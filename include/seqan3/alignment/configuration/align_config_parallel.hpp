// SPDX-FileCopyrightText: 2006-2025 Knut Reinert & Freie Universität Berlin
// SPDX-FileCopyrightText: 2016-2025 Knut Reinert & MPI für molekulare Genetik
// SPDX-License-Identifier: BSD-3-Clause

/*!\file
 * \brief Provides seqan3::align_cfg::parallel configuration.
 * \author Lydia Buntrock <lydia.buntrock AT fu-berlin.de>
 * \author Jörg Winkler <j.winkler AT fu-berlin.de>
 * \author Svenja Mehringer <svenja.mehringer AT fu-berlin.de>
 */

#pragma once

#include <seqan3/alignment/configuration/detail.hpp>
#include <seqan3/core/configuration/detail/configuration_element_parallel_mode.hpp>

namespace seqan3::align_cfg
{
/*!\brief Enables the parallel execution of the alignment algorithm if possible for the given configuration.
 * \ingroup alignment_configuration
 *
 * \details
 *
 * \include{doc} doc/fragments/alignment_configuration_align_config_parallel.md
 *
 * The value represents the number of threads to be used and must be greater than `0`.
 *
 * ### Example
 *
 * \include test/snippet/alignment/configuration/align_cfg_parallel_example.cpp
 *
 * \remark For a complete overview, take a look at \ref alignment_pairwise.
 */
using parallel = seqan3::detail::parallel_mode<
    std::integral_constant<seqan3::detail::align_config_id, seqan3::detail::align_config_id::parallel>>;

} // namespace seqan3::align_cfg
