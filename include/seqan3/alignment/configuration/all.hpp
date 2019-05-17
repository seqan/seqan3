// -----------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2019, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2019, Knut Reinert & MPI für molekulare Genetik
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at: https://github.com/seqan/seqan3/blob/master/LICENSE.md
// -----------------------------------------------------------------------------------------------------

 /*!\file
  * \brief Meta-header for the \link configuration alignment configuration module \endlink.
  * \author Rene Rahn <rene.rahn AT fu-berlin.de>
  */

 #pragma once

#include <seqan3/alignment/configuration/align_config_aligned_ends.hpp>
#include <seqan3/alignment/configuration/align_config_band.hpp>
#include <seqan3/alignment/configuration/align_config_edit.hpp>
#include <seqan3/alignment/configuration/align_config_gap.hpp>
#include <seqan3/alignment/configuration/align_config_max_error.hpp>
#include <seqan3/alignment/configuration/align_config_mode.hpp>
#include <seqan3/alignment/configuration/align_config_result.hpp>
#include <seqan3/alignment/configuration/align_config_scoring.hpp>
#include <seqan3/alignment/configuration/detail.hpp>

/*!\namespace seqan3::align_cfg
 * \brief A special sub namespace for the alignment configurations.
 */

/*!\defgroup alignment_configuration Configuration
 * \ingroup alignment
 *
 * ### Introduction
 *
 * In SeqAn the alignment algorithm can be configured in many different ways. The core of this configuration are the
 * different configuration elements that select specific features of the algorithm. To allow a maximal flexibility
 * the configuration is separated from the alignment interface. This means that before the alignment algorithm
 * (seqan3::align_pairwise) is invoked, the algorithm must be configured. The respective alignment configurations
 * live in their own namespace called seqan3::align_cfg. This namespace is used to disambiguate configurations for the
 * alignment algorithm with configurations from other algorithms.
 * To compute a pairwise alignment at least two configuration elements must be provided, namely the
 * the seqan3::align_cfg::mode and the seqan3::align_cfg::scoring.
 * The following code snippet shows the call of a pairwise alignment with the minimal configuration. It computes
 * a global alignment with custom scores for a match and a mismatch. The resulting score is `7`.
 *
 * \include test/snippet/alignment/configuration/minimal_alignment_config.cpp
 *
 * ### Combining configuration elements
 *
 * Configurations can be combined using the `|`-operator. If a combination is invalid, a static assertion is triggered
 * during compilation and will inform the user that the the last config cannot be combined with any of the configs from
 * the left-hand side of the configuration specification. Unfortunately, the names of the invalid
 * types cannot be printed within the static assert, but the following table shows which combinations are possible.
 * In general, the same configuration element cannot occur more than once inside of a configuration specification.
 *
 *<table class="wikitable" style="background-color:#ededed;font-size:85%;text-align:center;border: 1px solid black;border-collapse: collapse">
 *<tr>
 *<th style="border: 1px solid black;"> </th>
 *<th style="border: 1px solid black; vertical-align: middle; text-align: center; width: 7%;"> 0 </th>
 *<th style="border: 1px solid black; vertical-align: middle; text-align: center; width: 7%;"> 1 </th>
 *<th style="border: 1px solid black; vertical-align: middle; text-align: center; width: 7%;"> 2 </th>
 *<th style="border: 1px solid black; vertical-align: middle; text-align: center; width: 7%;"> 3 </th>
 *<th style="border: 1px solid black; vertical-align: middle; text-align: center; width: 7%;"> 4 </th>
 *<th style="border: 1px solid black; vertical-align: middle; text-align: center; width: 7%;"> 5 </th>
 *<th style="border: 1px solid black; vertical-align: middle; text-align: center; width: 7%;"> 6 </th>
 *<th style="border: 1px solid black; vertical-align: middle; text-align: center; width: 7%;"> 7 </th>
 *</tr>
 *<tr>
 *<th style="border: 1px solid black; text-align: left"> 0: seqan3::align_cfg::aligned_ends </th>
 *<td style="background: #ff9090; color: black; vertical-align: middle; text-align: center; border: 1px solid black;" class="table-no"> ✘ </td>
 *<td style="background: #90ff90; color: black; vertical-align: middle; text-align: center; border: 1px solid black;" class="table-yes"> ✓ </td>
 *<td style="background: #90ff90; color: black; vertical-align: middle; text-align: center; border: 1px solid black;" class="table-yes"> ✓ </td>
 *<td style="background: #90ff90; color: black; vertical-align: middle; text-align: center; border: 1px solid black;" class="table-yes"> ✓ </td>
 *<td style="background: #ff9090; color: black; vertical-align: middle; text-align: center; border: 1px solid black;" class="table-no"> ✘ </td>
 *<td style="background: #90ff90; color: black; vertical-align: middle; text-align: center; border: 1px solid black;" class="table-yes"> ✓ </td>
 *<td style="background: #90ff90; color: black; vertical-align: middle; text-align: center; border: 1px solid black;" class="table-yes"> ✓ </td>
 *<td style="background: #90ff90; color: black; vertical-align: middle; text-align: center; border: 1px solid black;" class="table-yes"> ✓ </td>
 *</tr>
 *<tr>
 *<th style="border: 1px solid black; text-align: left"> 1: seqan3::align_cfg::band </th>
 *<td style="background: #90ff90; color: black; vertical-align: middle; text-align: center; border: 1px solid black;" class="table-yes"> ✓ </td>
 *<td style="background: #ff9090; color: black; vertical-align: middle; text-align: center; border: 1px solid black;" class="table-no"> ✘ </td>
 *<td style="background: #90ff90; color: black; vertical-align: middle; text-align: center; border: 1px solid black;" class="table-yes"> ✓ </td>
 *<td style="background: #90ff90; color: black; vertical-align: middle; text-align: center; border: 1px solid black;" class="table-yes"> ✓ </td>
 *<td style="background: #90ff90; color: black; vertical-align: middle; text-align: center; border: 1px solid black;" class="table-yes"> ✓ </td>
 *<td style="background: #90ff90; color: black; vertical-align: middle; text-align: center; border: 1px solid black;" class="table-yes"> ✓ </td>
 *<td style="background: #90ff90; color: black; vertical-align: middle; text-align: center; border: 1px solid black;" class="table-yes"> ✓ </td>
 *<td style="background: #90ff90; color: black; vertical-align: middle; text-align: center; border: 1px solid black;" class="table-yes"> ✓ </td>
 *</tr>
 *<tr>
 *<th style="border: 1px solid black; text-align: left"> 2: seqan3::align_cfg::gap </th>
 *<td style="background: #90ff90; color: black; vertical-align: middle; text-align: center; border: 1px solid black;" class="table-yes"> ✓ </td>
 *<td style="background: #90ff90; color: black; vertical-align: middle; text-align: center; border: 1px solid black;" class="table-yes"> ✓ </td>
 *<td style="background: #ff9090; color: black; vertical-align: middle; text-align: center; border: 1px solid black;" class="table-no"> ✘ </td>
 *<td style="background: #90ff90; color: black; vertical-align: middle; text-align: center; border: 1px solid black;" class="table-yes"> ✓ </td>
 *<td style="background: #90ff90; color: black; vertical-align: middle; text-align: center; border: 1px solid black;" class="table-yes"> ✓ </td>
 *<td style="background: #90ff90; color: black; vertical-align: middle; text-align: center; border: 1px solid black;" class="table-yes"> ✓ </td>
 *<td style="background: #90ff90; color: black; vertical-align: middle; text-align: center; border: 1px solid black;" class="table-yes"> ✓ </td>
 *<td style="background: #90ff90; color: black; vertical-align: middle; text-align: center; border: 1px solid black;" class="table-yes"> ✓ </td>
 *</tr>
 *<tr>
 *<th style="border: 1px solid black; text-align: left"> 3: seqan3::global_alignment </th>
 *<td style="background: #90ff90; color: black; vertical-align: middle; text-align: center; border: 1px solid black;" class="table-yes"> ✓ </td>
 *<td style="background: #90ff90; color: black; vertical-align: middle; text-align: center; border: 1px solid black;" class="table-yes"> ✓ </td>
 *<td style="background: #90ff90; color: black; vertical-align: middle; text-align: center; border: 1px solid black;" class="table-yes"> ✓ </td>
 *<td style="background: #ff9090; color: black; vertical-align: middle; text-align: center; border: 1px solid black;" class="table-no"> ✘ </td>
 *<td style="background: #ff9090; color: black; vertical-align: middle; text-align: center; border: 1px solid black;" class="table-no"> ✘ </td>
 *<td style="background: #90ff90; color: black; vertical-align: middle; text-align: center; border: 1px solid black;" class="table-yes"> ✓ </td>
 *<td style="background: #90ff90; color: black; vertical-align: middle; text-align: center; border: 1px solid black;" class="table-yes"> ✓ </td>
 *<td style="background: #90ff90; color: black; vertical-align: middle; text-align: center; border: 1px solid black;" class="table-yes"> ✓ </td>
 *</tr>
 *<tr>
 *<th style="border: 1px solid black; text-align: left"> 4: seqan3::local_alignment </th>
 *<td style="background: #ff9090; color: black; vertical-align: middle; text-align: center; border: 1px solid black;" class="table-no"> ✘ </td>
 *<td style="background: #90ff90; color: black; vertical-align: middle; text-align: center; border: 1px solid black;" class="table-yes"> ✓ </td>
 *<td style="background: #90ff90; color: black; vertical-align: middle; text-align: center; border: 1px solid black;" class="table-yes"> ✓ </td>
 *<td style="background: #ff9090; color: black; vertical-align: middle; text-align: center; border: 1px solid black;" class="table-no"> ✘ </td>
 *<td style="background: #ff9090; color: black; vertical-align: middle; text-align: center; border: 1px solid black;" class="table-no"> ✘ </td>
 *<td style="background: #ff9090; color: black; vertical-align: middle; text-align: center; border: 1px solid black;" class="table-no"> ✘ </td>
 *<td style="background: #90ff90; color: black; vertical-align: middle; text-align: center; border: 1px solid black;" class="table-yes"> ✓ </td>
 *<td style="background: #90ff90; color: black; vertical-align: middle; text-align: center; border: 1px solid black;" class="table-yes"> ✓ </td>
 *</tr>
 *<tr>
 *<th style="border: 1px solid black; text-align: left"> 5: seqan3::align_cfg::max_error </th>
 *<td style="background: #90ff90; color: black; vertical-align: middle; text-align: center; border: 1px solid black;" class="table-yes"> ✓ </td>
 *<td style="background: #90ff90; color: black; vertical-align: middle; text-align: center; border: 1px solid black;" class="table-yes"> ✓ </td>
 *<td style="background: #90ff90; color: black; vertical-align: middle; text-align: center; border: 1px solid black;" class="table-yes"> ✓ </td>
 *<td style="background: #90ff90; color: black; vertical-align: middle; text-align: center; border: 1px solid black;" class="table-yes"> ✓ </td>
 *<td style="background: #ff9090; color: black; vertical-align: middle; text-align: center; border: 1px solid black;" class="table-no"> ✘ </td>
 *<td style="background: #ff9090; color: black; vertical-align: middle; text-align: center; border: 1px solid black;" class="table-no"> ✘ </td>
 *<td style="background: #90ff90; color: black; vertical-align: middle; text-align: center; border: 1px solid black;" class="table-yes"> ✓ </td>
 *<td style="background: #90ff90; color: black; vertical-align: middle; text-align: center; border: 1px solid black;" class="table-yes"> ✓ </td>
 *</tr>
 *<tr>
 *<th style="border: 1px solid black; text-align: left"> 6: seqan3::align_cfg::result </th>
 *<td style="background: #90ff90; color: black; vertical-align: middle; text-align: center; border: 1px solid black;" class="table-yes"> ✓ </td>
 *<td style="background: #90ff90; color: black; vertical-align: middle; text-align: center; border: 1px solid black;" class="table-yes"> ✓ </td>
 *<td style="background: #90ff90; color: black; vertical-align: middle; text-align: center; border: 1px solid black;" class="table-yes"> ✓ </td>
 *<td style="background: #90ff90; color: black; vertical-align: middle; text-align: center; border: 1px solid black;" class="table-yes"> ✓ </td>
 *<td style="background: #90ff90; color: black; vertical-align: middle; text-align: center; border: 1px solid black;" class="table-yes"> ✓ </td>
 *<td style="background: #90ff90; color: black; vertical-align: middle; text-align: center; border: 1px solid black;" class="table-yes"> ✓ </td>
 *<td style="background: #ff9090; color: black; vertical-align: middle; text-align: center; border: 1px solid black;" class="table-no"> ✘ </td>
 *<td style="background: #90ff90; color: black; vertical-align: middle; text-align: center; border: 1px solid black;" class="table-yes"> ✓ </td>
 *</tr>
 *<tr>
 *<th style="border: 1px solid black; text-align: left"> 7: seqan3::align_cfg::scoring </th>
 *<td style="background: #90ff90; color: black; vertical-align: middle; text-align: center; border: 1px solid black;" class="table-yes"> ✓ </td>
 *<td style="background: #90ff90; color: black; vertical-align: middle; text-align: center; border: 1px solid black;" class="table-yes"> ✓ </td>
 *<td style="background: #90ff90; color: black; vertical-align: middle; text-align: center; border: 1px solid black;" class="table-yes"> ✓ </td>
 *<td style="background: #90ff90; color: black; vertical-align: middle; text-align: center; border: 1px solid black;" class="table-yes"> ✓ </td>
 *<td style="background: #90ff90; color: black; vertical-align: middle; text-align: center; border: 1px solid black;" class="table-yes"> ✓ </td>
 *<td style="background: #90ff90; color: black; vertical-align: middle; text-align: center; border: 1px solid black;" class="table-yes"> ✓ </td>
 *<td style="background: #90ff90; color: black; vertical-align: middle; text-align: center; border: 1px solid black;" class="table-yes"> ✓ </td>
 *<td style="background: #ff9090; color: black; vertical-align: middle; text-align: center; border: 1px solid black;" class="table-no"> ✘ </td>
 *</tr>
 *</table>
 */
