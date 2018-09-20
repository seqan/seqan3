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
 * \author Marcel Ehrhardt <marcel.ehrhardt AT fu-berlin.de>
 * \brief Meta-header for the composition submodule; includes all headers from alphabet/composition/.
 */

#include <seqan3/alphabet/composition/cartesian_composition.hpp>
#include <seqan3/alphabet/composition/union_composition.hpp>

/*!\defgroup composition Composition
 * \brief Provides data structures joining multiple alphabets into a single alphabet.
 * \ingroup alphabet
 *
 * \par Introduction
 *
 * Composition alphabets are special alphabets that allow you to combine existing alphabets into new ones. For example,
 * you can add new characters to existing alphabets by using seqan3::union_composition or combine alphabets with quality
 * information by using seqan3::cartesian_composition.
 *
 * We have currently two major composition alphabets:
 * * seqan3::cartesian_composition which roughly corresponds to the Cartesian product of the given types. It
 *   behaves similar to std::tuple, but it is specialised for alphabets.
 * * seqan3::union_composition which roughly corresponds to the Union of the given types. It behaves similar to
 *   std::variant, but it is specialised for alphabets.
 */
