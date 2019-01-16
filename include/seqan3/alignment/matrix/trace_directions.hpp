// -----------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2019, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2019, Knut Reinert & MPI für molekulare Genetik
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at: https://github.com/seqan/seqan3/blob/master/LICENSE
// -----------------------------------------------------------------------------------------------------

/*!\file
 * \author Marcel Ehrhardt <marcel.ehrhardt AT fu-berlin.de>
 * \brief Contains the declaration of seqan3::detail::trace_directions.
 */

#pragma once

#include <seqan3/core/add_enum_bitwise_operators.hpp>

namespace seqan3::detail
{

/*!\brief The possible directions a trace can have. The values can be combined by the logical `|`-operator.
 * \ingroup alignment_matrix
 * \sa seqan3::add_enum_bitwise_operators <seqan3::detail::trace_directions> enables combining enum values.
 * \sa seqan3::detail::alignment_trace_matrix implementations use this enum as matrix entry type.
 */
enum struct trace_directions : uint8_t
{
    //!\brief No trace
    none      = 0b0000,
    //!\brief Trace comes from the diagonal entry.
    diagonal  = 0b0001,
    //!\brief Trace comes from the above entry.
    up        = 0b0010,
    //!\brief Trace comes from the left entry.
    left      = 0b0100
};

} // namespace seqan3::detail

namespace seqan3
{
//!\brief Enable bitwise operators for enum seqan3::detail::trace_directions.
//!\ingroup alignment_matrix
template <>
constexpr bool add_enum_bitwise_operators<seqan3::detail::trace_directions> = true;
} // namespace seqan3
