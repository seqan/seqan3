// -----------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2019, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2019, Knut Reinert & MPI für molekulare Genetik
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at: https://github.com/seqan/seqan3/blob/master/LICENSE
// -----------------------------------------------------------------------------------------------------

 /*!\file
  * \brief Meta-header for the \link configuration alignment configuration module \endlink.
  * \author Rene Rahn <rene.rahn AT fu-berlin.de>
  */

 #pragma once

/*!\defgroup configuration Configuration
 * \brief Data structures and utility functions for configuring alignment algorithm.
 * \ingroup alignment
 *
 * \todo Write detailed landing page.
 */
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
