// -----------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2019, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2019, Knut Reinert & MPI für molekulare Genetik
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at: https://github.com/seqan/seqan3/blob/master/LICENSE
// -----------------------------------------------------------------------------------------------------

/*!\file
 * \brief Provides seqan3::detail::affine_matrix_cell.
 * \author Rene Rahn <rene.rahn AT fu-berlin.de>
 */

#pragma once

#include <seqan3/alignment/matrix/trace_directions.hpp>
#include <seqan3/core/concept/core_language.hpp>
#include <seqan3/core/simd/concept.hpp>
#include <seqan3/std/concepts>

namespace seqan3::detail
{

/*!\brief A specific matrix entry for the affine alignment matrix.
 * \tparam score_type The type of the scoring matrix; must model seqan3::Arithmetic or seqan3::simd_concept.
 * \tparam trace_type The type of the trace matrix; must be of type seqan3::detail::trace_directions or model
 *                    seqan3::simd_concept and its element type is seqan3::detail::trace_directions or is
 *                    seqan3::detail::ignore_t to disable trace back computation.
 */
template <typename score_type, typename trace_type>
//!\cond
    requires (Arithmetic<score_type> || simd_concept<score_type>) &&
             (std::Same<trace_type, trace_directions> ||  simd_concept<trace_type> || decays_to_ignore_v<trace_type>)
//!\endcond
struct affine_matrix_cell
{
    score_type first; //!\< The first element.
    score_type second; //!\< The second element.
    trace_type dir; //!\< The trace direction.
};

/*!\name Type deduction guide
 * \{
 */

//!\brief Deduces the arguments from the constructor argument.
template <typename score_type, typename trace_type>
affine_matrix_cell(score_type, score_type, trace_type) -> affine_matrix_cell<score_type, trace_type>;
//!\}

} // namespace seqan3::detail
