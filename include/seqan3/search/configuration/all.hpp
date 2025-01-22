// SPDX-FileCopyrightText: 2006-2025 Knut Reinert & Freie Universität Berlin
// SPDX-FileCopyrightText: 2016-2025 Knut Reinert & MPI für molekulare Genetik
// SPDX-License-Identifier: BSD-3-Clause

/*!\file
 * \brief Meta-header for the \link search_configuration search configuration module \endlink.
 * \author Christopher Pockrandt <christopher.pockrandt AT fu-berlin.de>
 * \author Lydia Buntrock <lydia.buntrock AT fu-berlin.de>
 */

/*!\namespace seqan3::search_cfg
 * \brief A special sub namespace for the search configurations.
 */

/*!\defgroup search_configuration Configuration
 * \ingroup search
 * \see search
 * \brief Data structures and utility functions for configuring search algorithm.
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
 * "strategy", and which information should the \ref search_configuration_subsection_output "result" contain.
 * These configurations exist in their own namespace, namely seqan3::search_cfg, to disambiguate them from the
 * configuration of other algorithms.
 *
 * If no configuration is provided upon invoking the seqan3::search algorithm, a default configuration is provided:
 * \include test/snippet/search/configuration_default.cpp
 *
 * \section search_configuration_section_overview Overview on search configurations
 *
 * Configurations can be combined using the `|`-operator. If a combination is invalid, a static assertion is raised
 * during the compilation of the program. It will inform the user that some configurations cannot be combined together
 * into one search configuration. In general, the same configuration element cannot occur more than once inside of
 * a configuration specification. The following table shows which combinations are possible.
 *
 * | **Configuration group**                                                     | **0** | **1** | **2** | **3** | **4** | **5** | **6** |
 * |:----------------------------------------------------------------------------|:-----:|:-----:|:-----:|:-----:|:-----:|:-----:|:-----:|
 * | \ref seqan3::search_cfg::max_error_total  "0: Max error total"              |  ❌   |   ✅   |  ✅   |  ✅   |   ✅   |  ✅   |  ✅   |
 * | \ref seqan3::search_cfg::max_error_substitution "1: Max error substitution" |  ✅   |   ❌   |  ✅   |  ✅   |   ✅   |  ✅   |  ✅   |
 * | \ref seqan3::search_cfg::max_error_insertion "2: Max error insertion"       |  ✅   |   ✅   |  ❌   |  ✅   |   ✅   |  ✅   |  ✅   |
 * | \ref seqan3::search_cfg::max_error_deletion "3: Max error deletion"         |  ✅   |   ✅   |  ✅   |  ❌   |   ✅   |  ✅   |  ✅   |
 * | \ref search_configuration_subsection_output "4: Output"                     |  ✅   |   ✅   |  ✅   |  ✅   |   ❌   |  ✅   |  ✅   |
 * | \ref search_configuration_subsection_hit_strategy "5: Hit"                  |  ✅   |   ✅   |  ✅   |  ✅   |   ✅   |  ❌   |  ✅   |
 * | \ref seqan3::search_cfg::parallel "6: Parallel"                             |  ✅   |   ✅   |  ✅   |  ✅   |   ✅   |  ✅   |  ❌   |
 *
 * \subsection search_configuration_subsection_error 0 - 3: Max Error Configuration
 *
 * This configuration can be used to specify the number or rate of error types.
 * It restricts the number of substitutions, insertions, deletions and total errors within the search to the given
 * values. A mismatch corresponds to diverging bases between text and query for a certain position. An insertion
 * corresponds to a base inserted into the query that does not occur in the text at the respective position. A deletion
 * corresponds to a base deleted from the query sequence that does occur in the indexed text. Deletions at the beginning
 * and at the end of the sequence are not considered during a search.
 *
 * The following rules apply when selecting the max error configuration:
 * First, if seqan3::search_cfg::max_error_total is specified, then all error types are set to
 * the value of the total error configuration. For any other specified error configuration the value is set accordingly,
 * but will not exceed the total error if given. For example, if a configuration profile sets the total max error to 3
 * and the insertion error to 1, then the search will at most consider one insertion but allow up to 3 deletions
 * and 3 substitutions during the search, while allowing at most 3 errors in total.
 * On the other hand, if the total error is not specified in the search configuration, it will be set to the sum of the
 * other configurations. This means that in the default case all errors are set to 0 and therefore an exact search is
 * conducted.
 *
 * The configuration elements can be initialised by an absolute error count or an error rate:
 * | seqan3::search_cfg::max_error_*¹| Behaviour                                                       |
 * |---------------------------------|-----------------------------------------------------------------|
 * | seqan3::search_cfg::error_rate  | Specify the error rate (\f$\in [0,1]\f$).                       |
 * | seqan3::search_cfg::error_count | Specify a discrete number of allowed errors (\f$\mathbb{W}\f$). |
 *
 * ¹: max_error_total, max_error_substitution, max_error_insertion, max_error_deletion
 *
 * ### Example
 *
 * \include test/snippet/search/configuration_error.cpp
 *
 * \subsection search_configuration_subsection_output 4. Output Configuration
 *
 * The seqan3::search interface returns a seqan3::algorithm_result_generator_range. This range is a lazy single
 * pass input range over the computed hits and the range's element types are seqan3::search_result objects.
 * Even if only a single query is searched, a range will be returned since it could be possible that
 * one search produces multiple hits, e.g. to find all best hits.
 * The following output configurations exists:
 *
 * * seqan3::search_cfg::output_query_id
 * * seqan3::search_cfg::output_reference_id
 * * seqan3::search_cfg::output_reference_begin_position
 * * seqan3::search_cfg::output_index_cursor
 *
 * Each of the configuration elements corresponds to a member access function inside of the returned
 * seqan3::search_result object. If you do not specify any output configuration in the search configuration, then the
 * default output contains the query and reference id as well as the reference begin position. If you customise the output
 * configuration, then only those are available in the final seqan3::search_result that are specified.
 *
 * \note If you try to call a function for an entity that was not configured a static
 *       assertion will be raised during compilation of the program.
 *
 * \include test/snippet/search/configuration_output.cpp
 *
 * The index cursor is an advanced data structure that lets you navigate within the index.
 * See seqan3::fm_index_cursor and seqan3::bi_fm_index_cursor for more information.
 * If you don't need the reference id nor the position, returning only the cursor is faster.
 * This is, because the operation to get the id and position of a hit can be computationally intensive
 * depending on the underlying index structure.
 *
 * \note A single index cursor points to **a range of text positions**. Although the normal use case is to return
 *       either the cursor or the positions, both can be returned simultaneously. In this case, the same cursor will
 *       be copied into the seqan3::search_result for each of its associated positions.
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
 * Instead, one can use the dynamic hit configuration element seqan3::search_cfg::hit.
 * This configuration element allows to set one of the above mentioned hit configurations at runtime. Later during the
 * configuration phase of the search algorithm the selected search configuration is used for the final search algorithm.
 * If the dynamic hit configuration is default constructed it does not hold any hit configuration. If you call search
 * with the dynamic configuration in this state an exception will be thrown.
 * Also note that using the dynamic configuration might have implications on the compile time, so we recommend to use
 * the static configurations if only a single hit strategy is supported.
 * The following example demonstrates the usage of the dynamic configuration:
 *
 * \include test/snippet/search/dynamic_hit_configuration_example.cpp
 *
 * \subsection search_configuration_subsection_parallel 6: Parallel Configuration
 *
 * This configuration determines the maximal number of threads the search algorithm can use.
 *
 * The seqan3::search_cfg::parallel configuration element can be combined with any other search configuration.
 *
 * \include test/snippet/search/configuration_parallel.cpp
 *
 * ### User callback
 *
 * In the default case, a call to seqan3::search returns a lazy range over the results of the search. This lazy range
 * has the advantage that the results are always in a deterministic order even if the search is executed in parallel.
 * Sometimes, however, it might be desirable to provide a user defined callback.
 * To do so, one can use the configuration element seqan3::search_cfg::on_result. This configuration element
 * is initialised with a user defined callback, e.g. a lambda function, which will be invoked with a generated
 * seqan3::search_result whenever a hit was found.
 * This has two implications. First, the return type of the seqan3::search function changes to `void`, i.e. it
 * returns nothing. Second, in a parallel execution of the search, the order of the hits is not deterministic and the
 * user has to make sure that concurrent invocations of the given callback are safe.
 *
 * The following snippet demonstrates the basic use case for this configuration element:
 *
 * \include test/snippet/search/search_with_user_callback.cpp
 */

#pragma once

#include <seqan3/search/configuration/default_configuration.hpp>
#include <seqan3/search/configuration/hit.hpp>
#include <seqan3/search/configuration/max_error.hpp>
#include <seqan3/search/configuration/max_error_common.hpp>
#include <seqan3/search/configuration/on_result.hpp>
#include <seqan3/search/configuration/output.hpp>
#include <seqan3/search/configuration/parallel.hpp>
#include <seqan3/search/configuration/result_type.hpp>
