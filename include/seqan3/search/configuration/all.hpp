// -----------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2020, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2020, Knut Reinert & MPI für molekulare Genetik
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at: https://github.com/seqan/seqan3/blob/master/LICENSE.md
// -----------------------------------------------------------------------------------------------------

/*!\file
 * \brief Meta-header for the \link search_configuration search configuration module \endlink.
 * \author Christopher Pockrandt <christopher.pockrandt AT fu-berlin.de>
 * \author Lydia Buntrock <lydia.buntrock AT fu-berlin.de>
 */

#pragma once

#include <seqan3/core/algorithm/configuration.hpp>
#include <seqan3/search/configuration/default_configuration.hpp>
#include <seqan3/search/configuration/detail.hpp>
#include <seqan3/search/configuration/max_error.hpp>
#include <seqan3/search/configuration/hit.hpp>
#include <seqan3/search/configuration/output.hpp>
#include <seqan3/search/configuration/parallel.hpp>
#include <seqan3/search/configuration/result_type.hpp>

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
 * \section search_configuration_section_introduction Introduction
 *
 * In SeqAn, the search algorithm uses a configuration object to determine the desired amount of
 * \ref seqan3::search_cfg::max_error_total "total errors",
 * of \ref seqan3::search_cfg::max_error_substitution "substitution errors",
 * of \ref seqan3::search_cfg::max_error_insertion "insertion errors",
 * and of \ref seqan3::search_cfg::max_error_deletion "deletion errors",
 * where all can be given as an \ref seqan3::search_cfg::error_count "absolute number"
 * or a \ref seqan3::search_cfg::error_rate "rate" of \ref search_configuration_subsection_error "errors".
 * Furthermore, it can be configured what hits are reported based on a \ref search_configuration_subsection_hit_strategy
 * "strategy", and how to \ref seqan3::search_cfg::output "output" the results.
 * These configurations exist in their own namespace, namely seqan3::search_cfg, to disambiguate them from the
 * configuration of other algorithms.
 *
 * If no configuration is provided upon invoking the seqan3::search algorithm, a default configuration is provided:
 * \include test/snippet/search/configuration_default.cpp
 *
 * \section search_configuration_section_overview Overview on search configurations
 *
 * Configurations can be combined using the `|`-operator. If a combination is invalid, a static assertion is triggered
 * during compilation and will inform the user that the the last config cannot be combined with any of the configs from
 * the left-hand side of the configuration specification. Unfortunately, the names of the invalid
 * types cannot be printed within the static assert, but the following table shows which combinations are possible.
 * In general, the same configuration element cannot occur more than once inside of a configuration specification.
 *
 * | **Configuration group**                                                     | **0** | **1** | **2** | **3** | **4** | **5** | **6** |
 * | ----------------------------------------------------------------------------|-------|-------|-------|-------|-------|-------|-------|
 * | \ref seqan3::search_cfg::max_error_total  "0: Max error total"              |  ❌   |  ✅   |  ✅   |  ✅   |  ✅   |  ✅   |  ✅   |
 * | \ref seqan3::search_cfg::max_error_substitution "1: Max error substitution" |  ✅   |  ❌   |  ✅   |  ✅   |  ✅   |  ✅   |  ✅   |
 * | \ref seqan3::search_cfg::max_error_insertion "2: Max error insertion"       |  ✅   |  ✅   |  ❌   |  ✅   |  ✅   |  ✅   |  ✅   |
 * | \ref seqan3::search_cfg::max_error_deletion "3: Max error deletion"         |  ✅   |  ✅   |  ✅   |  ❌   |  ✅   |  ✅   |  ✅   |
 * | \ref seqan3::search_cfg::output "4: Output"                                 |  ✅   |  ✅   |  ✅   |  ✅   |  ❌   |  ✅   |  ✅   |
 * | \ref search_configuration_subsection_hit_strategy "5: Hit"                  |  ✅   |  ✅   |  ✅   |  ✅   |  ✅   |  ❌   |  ✅   |
 * | \ref seqan3::search_cfg::parallel "6: Parallel"                             |  ✅   |  ✅   |  ✅   |  ✅   |  ✅   |  ✅   |  ❌   |
 *
 * \subsection search_configuration_search_result Search result type
 *
 * \copydetails seqan3::search_result
 *
 * \subsection search_configuration_subsection_error 0 - 3: Max Error Configuration
 *
 * This configuration can be used to specify the number or rate of error types.
 * It restricts the number of substitutions, insertions, deletions and total errors within the search to the given
 * values. A mismatch corresponds to diverging bases between text and query for a certain position. An insertion
 * corresponds to a base inserted into the query that does not occur in the text at the position. A deletion corresponds
 * to a base deleted from the query sequence that does occur in the indexed text. Deletions at the beginning and at the
 * end of the sequence are not considered during a search.
 *
 * They will behave as follows:
 * If seqan3::search_cfg::total and any other error type are specified, all types are set to the respective values.
 * If one or more other error types are configured, but no total, then total is set to the sum of the error types.
 *
 * These configuration elements can be given by a number or rate of errors:
 * | seqan3::search_cfg::max_error_*¹| Behaviour                                    |
 * |---------------------------------|----------------------------------------------|
 * | seqan3::search_cfg::error_rate  | Specify the error rate.                      |
 * | seqan3::search_cfg::error_count | Specify a descrete number of allowed errors. |
 *
 * ¹: max_error_total, max_error_substitution, max_error_insertion, max_error_deletion
 *
 * ### Example
 *
 * \include test/snippet/search/configuration_error.cpp
 *
 * \subsection search_configuration_subsection_hit_strategy 5: Hit Configuration
 *
 * This configuration can be used to determine which hits are reported.
 * Currently these strategies are available:
 *
 * | Hit Configurations                  | Behaviour                                                           |
 * |-------------------------------------|---------------------------------------------------------------------|
 * | seqan3::search_cfg::hit_all         | Report all hits within error bounds.                                |
 * | seqan3::search_cfg::hit_all_best    | Report all hits with the lowest number of errors within the bounds. |
 * | seqan3::search_cfg::hit_single_best | Report one best hit (hit with lowest error) within bounds.          |
 * | seqan3::search_cfg::hit_strata      | Report all hits within best + `stratum` errors.                     |
 *
 * The individual configuration elements to select a search strategy cannot be combined with each other
 * (mutual exclusivity).
 *
 * \include test/snippet/search/hit_configuration_examples.cpp
 *
 * ### Dynamic hit configuration
 *
 * Sometimes a program needs to support different hit strategies based on some user input. Since these are mostly
 * runtime decisisons the code can become quite cumbersome to handle the static hit configurations.
 * Instead, one can use the dynamic hit configuration element seqan3::search::cfg::hit.
 * This configuration element allows to set one of the above mentioned hit configurations at runtime. Later during the
 * configuration phase of the search algorithm the selected search configuration is used for the final search algorithm.
 * If the dynamic hit configuration is default constructed it does not hold any hit configuration. If you call search
 * with the dynamic configuration in this state an exception will be thrown.
 * Also note that using the dynamic configuration might have implications on the compile time, so we recommend to use
 * the static configurations if only a single hit strategy is supported.
 * The following example demonstrates the usage of the dynamic configuration:
 *
 * \include test/snippet/search/dynamic_hit_configuration_example.cpp
 */
