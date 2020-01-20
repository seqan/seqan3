// -----------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2020, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2020, Knut Reinert & MPI für molekulare Genetik
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at: https://github.com/seqan/seqan3/blob/master/LICENSE.md
// -----------------------------------------------------------------------------------------------------

/*!\file
 * \brief Provides seqan3::detail::alignment_trace_matrix_base.
 * \author Rene Rahn <rene.rahn AT fu-berlin.de>
 */

#pragma once

#include <vector>

#include <seqan3/alignment/matrix/alignment_coordinate.hpp>
#include <seqan3/alignment/matrix/detail/two_dimensional_matrix.hpp>
#include <seqan3/alignment/matrix/trace_directions.hpp>
#include <seqan3/core/simd/concept.hpp>
#include <seqan3/range/container/aligned_allocator.hpp>

namespace seqan3::detail
{

/*!\brief A crtp-base class for alignment traceback matrices.
 * \tparam derived_t The derived type.
 * \tparam trace_t   The type of the trace directions.
 * \ingroup alignment_matrix
 *
 * \details
 *
 * Manages the actual storage as a std::vector. How much memory is allocated is handled by the derived type.
 * The `trace_t` must be either a seqan3::detail::trace_directions enum value or a seqan3::detail::simd_conceptvector
 * over seqan3::detail::trace_directions.
 */
template <typename trace_t>
struct alignment_trace_matrix_base
{
protected:
    static_assert(std::same_as<trace_t, trace_directions> || simd_concept<trace_t>,
                  "Value type must either be a trace_directions object or a simd vector.");

    //!\brief The coordinate type.
    using coordinate_type = advanceable_alignment_coordinate<advanceable_alignment_coordinate_state::row>;
    //!\brief The actual element type.
    using element_type = trace_t;
    //!\brief The allocator type. Uses seqan3::aligned_allocator if storing seqan3::detail::simd_concepttypes.
    using allocator_type = std::conditional_t<detail::simd_concept<trace_t>,
                                              aligned_allocator<element_type, sizeof(element_type)>,
                                              std::allocator<element_type>>;
    //!\brief The type of the underlying memory pool.
    using pool_type = two_dimensional_matrix<element_type, allocator_type, matrix_major_order::column>;
    //!\brief The size type.
    using size_type = size_t;

public:
    //!\brief The linearised matrix storing the trace data in column-major-order.
    pool_type data{};
    //!\brief Internal cache for the trace values to the left.
    std::vector<element_type, allocator_type> cache_left{};
    //!\brief Internal cache for the last trace value above.
    element_type cache_up{};
    //!\brief The number of columns.
    size_type num_cols{};
    //!\brief The number of num_rows.
    size_type num_rows{};
};

} // namespace seqan3::detail
