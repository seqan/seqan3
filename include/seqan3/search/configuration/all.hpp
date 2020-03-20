// -----------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2020, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2020, Knut Reinert & MPI für molekulare Genetik
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at: https://github.com/seqan/seqan3/blob/master/LICENSE.md
// -----------------------------------------------------------------------------------------------------

/*!\file
 * \brief Meta-header for the \link search_configuration search configuration module \endlink.
 * \author Christopher Pockrandt <christopher.pockrandt AT fu-berlin.de>
 */

#pragma once

#include <seqan3/core/algorithm/configuration.hpp>
#include <seqan3/search/configuration/default_configuration.hpp>
#include <seqan3/search/configuration/detail.hpp>
#include <seqan3/search/configuration/max_error.hpp>
#include <seqan3/search/configuration/max_error_rate.hpp>
#include <seqan3/search/configuration/mode.hpp>
#include <seqan3/search/configuration/output.hpp>
#include <seqan3/search/configuration/parallel.hpp>

/*!\namespace seqan3::search_cfg
 * \brief A special sub namespace for the search configurations.
 */

/*!\defgroup search_configuration Configuration
 * \ingroup search
 * \brief Data structures and utility functions for configuring search algorithm.
 * \see search
 *
 * \details
 *
 * ### Introduction
 *
 * In SeqAn the search algorithm uses a configuration object to determine the desired
 * \ref seqan3::search_cfg::max_error "number"/\ref seqan3::search_cfg::max_error_rate "rate" of errors,
 * what hits are considered as \ref seqan3::search_cfg::mode "results", and how to
 * \ref seqan3::search_cfg::output "output" the result.
 * These configurations exist in their own namespace, namely seqan3::search_cfg, to disambiguate them from the
 * configuration of other algorithms.
 *
 * If no configuration is provided upon invoking the seqan3::search algorithm, a default configuration is provided:
 * \include test/snippet/search/configuration_default.cpp
 *
 * ### Combining configuration elements
 *
 * Configurations can be combined using the `|`-operator. If a combination is invalid, a static assertion is triggered
 * during compilation and will inform the user that the the last config cannot be combined with any of the configs from
 * the left-hand side of the configuration specification. Unfortunately, the names of the invalid
 * types cannot be printed within the static assert, but the following table shows which combinations are possible.
 * In general, the same configuration element cannot occur more than once inside of a configuration specification.
 *
 * | **Config**                                                  | **0** | **1** | **2** | **3** | **4** |
 * | ------------------------------------------------------------|-------|-------|-------|-------|-------|
 * | \ref seqan3::search_cfg::max_error  "0: Max error"          |  ❌   |  ❌   |  ✅   |  ✅   |  ✅   |
 * | \ref seqan3::search_cfg::max_error_rate "1: Max error rate" |  ❌   |  ❌   |  ✅   |  ✅   |  ✅   |
 * | \ref seqan3::search_cfg::output "2: Output"                 |  ✅   |  ✅   |  ❌   |  ✅   |  ✅   |
 * | \ref seqan3::search_cfg::mode "3: Mode"                     |  ✅   |  ✅   |  ✅   |  ❌   |  ✅   |
 * | \ref seqan3::search_cfg::parallel "4: Parallel"             |  ✅   |  ✅   |  ✅   |  ✅   |  ❌   |
 */
