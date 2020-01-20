// -----------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2020, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2020, Knut Reinert & MPI für molekulare Genetik
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at: https://github.com/seqan/seqan3/blob/master/LICENSE.md
// -----------------------------------------------------------------------------------------------------

//!\cond DEV
/*!\file
 * \brief Meta-header for the \link alignment_policy alignment policy submodule \endlink.
 * \author Rene Rahn <rene.rahn AT fu-berlin.de>
 */
//!\endcond

#pragma once

#include <seqan3/alignment/pairwise/policy/affine_gap_init_policy.hpp>
#include <seqan3/alignment/pairwise/policy/affine_gap_policy.hpp>
#include <seqan3/alignment/pairwise/policy/alignment_matrix_policy.hpp>
#include <seqan3/alignment/pairwise/policy/find_optimum_policy.hpp>
#include <seqan3/alignment/pairwise/policy/scoring_scheme_policy.hpp>
#include <seqan3/alignment/pairwise/policy/simd_affine_gap_policy.hpp>
#include <seqan3/alignment/pairwise/policy/simd_find_optimum_policy.hpp>

//!\cond DEV
/*!\defgroup alignment_policy Alignment policies
 * \ingroup pairwise_alignment
 * \brief Provides policies for the alignment algorithm.
 *
 * # Introduction
 *
 * The standard pairwise alignment algorithm in SeqAn is implemented in many variations. It supports the standard
 * algorithms such as the global, local, and semi-global alignment using different kinds of scoring matrices for
 * nucleotide or amino acid alphabets. It further allows for computing only the score or the begin and end positions
 * or even the traceback. In addition, the algorithms can be executed in highly parallel environments using SIMD
 * (Single Instruction Multiple Data) vectorisation and multi-threading. The combination of all of these variations
 * leads to a huge number of different implementations of the same algorithm. Hence it is desirable to reduce the
 * code duplication in order to increase maintenance and extension of the alignment algorithms in the future. To
 * achieve this the alignment algorithm type is specialised with alignment policies.
 *
 * ### Policy state
 *
 * Policies can have an internal state to manage some additional variables. However, be careful with using
 * non-stateless policies, as they can affect the internal memory layout which can result in performance regressions.
 * The state of a policy should therefore be carefully tested with benchmarks.
 *
 * ### Customising the alignment algorithm
 *
 * An alignment policy serves as a customisation point to the alignment algorithm which has to implement a specific
 * set of functions that are called by the actual seqan3::detail::alignment_algorithm type. These policies further
 * separate logical units of the alignment algorithm, i.e. the initialisation, the computation, and the memory
 * allocation of the alignment matrix.
 *
 * # Gap policies
 *
 * Gap policies are used to initialise and to compute the cells within the alignment matrix.
 * The gap policies are further divided into a policy initialising the matrix and computing the cells.
 * The following table shows the functions that need to be implemented by a gap policy.
 *
 * Function name     | Arguments                              | Return value                         |
 * ----------------- | -------------------------------------- | ------------------------------------ |
 * `compute_cell`    | `cell &&`, `cache &`, `score const &`  | `void`                               |
 * `make_cache`      | `gap_scheme const &`                   | alignment_algorithm_state            |
 *
 * - compute_cell:
 *
 *    This function implements the kernel that computes the score and, if enabled, the traceback direction.
 *    It gets as an input the dereferenced value of the scoring matrix (see above), the current alignment algorithm
 *    cache and the score of comparing two letter of the used alphabet using the passed scoring scheme.
 *
 * - make_cache:
 *
 *     Initialises and returns the cache used in the alignment algorithm. It must return a type that allows structured
 *     bindings, with the first argument being of the same type as the allocator value type used in the matrix policy,
 *     the second argument being the costs for gap opening + gap and the third argument being the costs for a gap.
 *
 * ### Existing gap policies:
 *
 *  - seqan3::detail::affine_gap_policy
 *
 * The following table displays requirements for the corresponding gap intialisation policy:
 *
 * Function name      | Arguments                        | Return value                         |
 * ------------------ | -------------------------------- | ------------------------------------ |
 * `init_origin_cell` | `cell &&`, `cache &`             | `void`                               |
 * `init_column_cell` | `cell &&`, `cache &`             | `void`                               |
 * `init_row_cell`    | `cell &&`, `cache &`             | `void`                               |
 *
 * - init_origin_cell:
 *
 *     Implements the initialisation of the matrix origin \f$ (M(0,0)) \f$. This function gets the current cell
 *     dereferenced from the score matrix and the alignment algorithm cache.
 *
 * - init_column_cell:
 *
 *     Implements the initialisation of the cells in the first row of the matrix \f$ (M(i,0)) \f$.
 *     This function gets the current cell dereferenced from the score matrix and the alignment algorithm cache.
 *
 * - init_row_cell:
 *
 *     Implements the initialisation of the cells in the first column of the matrix \f$ (M(0,j)) \f$.
 *     This function gets the current cell dereferenced from the score matrix and the alignment algorithm cache.
 *
 * ### Existing gap init policies:
 *
 *  - seqan3::detail::affine_gap_init_policy
 *
 * # Find optimum policies
 *
 * These policies are used to define the search space of the alignment optimum.
 *
 * Function name             | Arguments                     | Return value                         |
 * ------------------------- | ----------------------------- | ------------------------------------ |
 * `check_score`             | `cell const`, `optimum &`     | `void`                               |
 * `check_score_last_row`    | `cell const`, `optimum &`     | `void`                               |
 * `check_score_last_column` | `cell const`, `optimum &`     | `void`                               |
 *
 * - check_score:
 *
 *    Is called for every cell in the dynamic programming matrix. Might be a "NO-OP".
 *    The left operand is the score of the current cell and the right operand is the current optimum that will
 *    be updated if the current score was higher.
 *
 * - check_score_last_row:
 *
 *    Is called only for the cells in the last row of the dynamic programming matrix. Might be a NO-OP.
 *    The left operand is the score of the current cell and the right operand is the current optimum that will
 *    be updated if the current score was higher.
 *
 * - check_score_last_column:
 *
 *    Is called only for the cells in the last column of the dynamic programming matrix.
 *    The left operand is the entire last column, because the algorithm's layout is column based.
 *    This function might skip the entire range and only compare the last value (the score for the global alignment)
 *    depending on the alignment configuration.
 *
 * ### Existing optimum policies:
 *
 *  - seqan3::detail::find_optimum_policy
 */
//!\endcond
