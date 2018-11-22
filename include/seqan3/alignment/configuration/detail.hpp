// ============================================================================
//                 SeqAn - The Library for Sequence Analysis
// ============================================================================
//
// Copyright (c) 2006-2018, Knut Reinert & Freie Universitaet Berlin
// Copyright (c) 2016-2018, Knut Reinert & MPI Molekulare Genetik
// All rights reserved.
//
// Redistribution and use in source and binary forms, with or without
// modification, are permitted provided that the following conditions are met:
//
//     * Redistributions of source code must retain the above copyright
//       notice, this list of conditions and the following disclaimer.
//     * Redistributions in binary form must reproduce the above copyright
//       notice, this list of conditions and the following disclaimer in the
//       documentation and/or other materials provided with the distribution.
//     * Neither the name of Knut Reinert or the FU Berlin nor the names of
//       its contributors may be used to endorse or promote products derived
//       from this software without specific prior written permission.
//
// THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
// AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
// IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
// ARE DISCLAIMED. IN NO EVENT SHALL KNUT REINERT OR THE FU BERLIN BE LIABLE
// FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL
// DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR
// SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER
// CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT
// LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY
// OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH
// DAMAGE.
//
// ============================================================================

/*!\file
 * \brief Provides some utility functions for the alignment configurations.
 * \author Rene Rahn <rene.rahn AT fu-berlin.de>
 */

#pragma once

#include <seqan3/core/algorithm/configuration_utility.hpp>

namespace seqan3::detail
{

/*!\brief An internal enum to check for a consistent configuration object.
 * \ingroup configuration
 */
enum struct align_config_id : uint8_t
{
    aligned_ends,
    band,
    gap,
    global,
    max_error,
    result,
    scoring,
    SIZE
};

// ----------------------------------------------------------------------------
// compatibility_table
// ----------------------------------------------------------------------------

/*!\brief Declaration of algorithm specific compatibility table.
 * \ingroup algorithm
 * \tparam algorithm_id_type The type of the algorithm specific id. Algorithm configurations must maintain
 *                           this table to allow validation checks.
 */
template <>
inline constexpr std::array<std::array<bool, static_cast<uint8_t>(align_config_id::SIZE)>,
                            static_cast<uint8_t>(align_config_id::SIZE)> compatibility_table<align_config_id>
{
{   //0  1  2  3  4  5  6
    { 0, 1, 1, 1, 1, 1, 1}, // 0: aligned_ends
    { 1, 0, 1, 1, 1, 1, 1}, // 1: band
    { 1, 1, 0, 1, 1, 1, 1}, // 2: gap
    { 1, 1, 1, 0, 1, 1, 1}, // 3: global
    { 1, 1, 1, 1, 0, 1, 1}, // 4: max_error
    { 1, 1, 1, 1, 1, 0, 1}, // 5: result
    { 1, 1, 1, 1, 1, 1, 0}  // 6: scoring
}
};

} // namespace seqan3::detail
