// SPDX-FileCopyrightText: 2006-2025 Knut Reinert & Freie Universität Berlin
// SPDX-FileCopyrightText: 2016-2025 Knut Reinert & MPI für molekulare Genetik
// SPDX-License-Identifier: BSD-3-Clause

/*!\file
 * \brief Meta-header for the \link alignment_configuration Alignment / Configuration submodule \endlink.
 * \author Rene Rahn <rene.rahn AT fu-berlin.de>
 */

/*!\namespace seqan3::align_cfg
 * \brief A special sub namespace for the alignment configurations.
 */

/*!\if DEV
  * \namespace seqan3::align_cfg::detail
  * \copydoc seqan3::detail
  * \endif
  */

/*!\defgroup alignment_configuration Configuration
 * \ingroup alignment
 * \brief Provides configuration elements for the pairwise alignment configuration.
 *
 * See the detailed \ref alignment_pairwise documentation for more details.
 *
 * \see alignment_pairwise
 */

#pragma once

#include <seqan3/alignment/configuration/align_config_band.hpp>
#include <seqan3/alignment/configuration/align_config_debug.hpp>
#include <seqan3/alignment/configuration/align_config_edit.hpp>
#include <seqan3/alignment/configuration/align_config_gap_cost_affine.hpp>
#include <seqan3/alignment/configuration/align_config_method.hpp>
#include <seqan3/alignment/configuration/align_config_min_score.hpp>
#include <seqan3/alignment/configuration/align_config_on_result.hpp>
#include <seqan3/alignment/configuration/align_config_output.hpp>
#include <seqan3/alignment/configuration/align_config_parallel.hpp>
#include <seqan3/alignment/configuration/align_config_result_type.hpp>
#include <seqan3/alignment/configuration/align_config_score_type.hpp>
#include <seqan3/alignment/configuration/align_config_scoring_scheme.hpp>
#include <seqan3/alignment/configuration/align_config_vectorised.hpp>
