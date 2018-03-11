// ============================================================================
//                 SeqAn - The Library for Sequence Analysis
// ============================================================================
//
// Copyright (c) 2006-2017, Knut Reinert & Freie Universitaet Berlin
// Copyright (c) 2016-2017, Knut Reinert & MPI Molekulare Genetik
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
 * \brief Meta-header for the \link action action submodule \endlink.
 * \author Hannes Hauswedell <hannes.hauswedell AT fu-berlin.de>
 */

#pragma once

/*!\defgroup action Action
 * \brief Actions are "eager range combinators" that modify ranges in-place.
 * \ingroup range
 * \sa https://ericniebler.github.io/range-v3/index.html#range-actions
 * \sa range/action.hpp
 *
 * SeqAn3 makes use of actions as defined in the
 * [Ranges Technical Specification](http://en.cppreference.com/w/cpp/experimental/ranges). Currently the
 * implementation is based on the [range-v3 library](https://github.com/ericniebler/range-v3) and all those actions
 * are available in the namespace ranges::action, see
 * [the overview](https://ericniebler.github.io/range-v3/index.html#range-actions) for more details.
 *
 * This submodule provides additional actions, specifically for operations on biological data and
 * sequence analysis.
 *
 * \attention
 * To prevent naming conflicts, all SeqAn actions are inside the namespace seqan3::action.
 *
 * \par Example
 *
 * ```cpp
 * dna4_vector vec{"ACGGTC"_dna4};
 * vec |= action::complement;                           // == "TGCCAG"
 * vec |= ranges::action::reverse;                      // == "GACCGT"
 *
 * // or in one line:
 * vec = "ACGGTC"_dna4;                                 // == "ACGGTC"
 * vec |= action::complement | ranges::action::reverse; // == "GACCGT"
 * ```
 */

/*!
 * \namespace seqan3::action
 * \brief The SeqAn3 namespace for actions.
 * \sa action
 *
 * Since actions often have name clashes with regular functions and views they are implemented in the sub
 * namespace `action`.
 *
 * See the \link action action submodule \endlink of the range module for more details.
 */
