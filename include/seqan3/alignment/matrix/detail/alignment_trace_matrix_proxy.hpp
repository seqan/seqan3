// -----------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2020, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2020, Knut Reinert & MPI für molekulare Genetik
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at: https://github.com/seqan/seqan3/blob/master/LICENSE.md
// -----------------------------------------------------------------------------------------------------

/*!\file
 * \brief Provides seqan3::detail::alignment_trace_matrix_proxy.
 * \author Rene Rahn <rene.rahn AT fu-berlin.de>
 */

#pragma once

#include <seqan3/alignment/matrix/trace_directions.hpp>
#include <seqan3/core/simd/concept.hpp>
#include <seqan3/std/concepts>

namespace seqan3::detail
{

/*!\brief A proxy type for a unified access to the traceback matrix during alignment computation.
 * \tparam trace_type The type wrapped by this proxy.
 * \ingroup alignment_matrix
 *
 * \details
 *
 * Provides named accessors to the respective values of the traceback matrix during the alignment computation.
 * The `trace_type` must be either a seqan3::detail::trace_directions value or a seqan3::detail::simd_conceptvector
 * over seqan3::detail::trace_directions.
 */
template <typename coordinate_type, typename trace_type>
struct alignment_trace_matrix_proxy
{
    static_assert(std::same_as<trace_type, trace_directions> || simd_concept<trace_type> ||
                  decays_to_ignore_v<trace_type>,
                  "Value type must either be a trace_directions object, a simd vector over such or std::ignore.");

    coordinate_type coordinate{}; //!< The current coordinate.
    trace_type & current; //!< Reference to the current element.
    trace_type & r_left; //!< Reference to the element to the left.
    trace_type & w_left; //!< Reference to the next element to the left.
    trace_type & up; //!< Reference to the element above.
};

} // namespace seqan3::detail
