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
 * \ingroup view
 * \brief Meta-header for the \link view view submodule \endlink.
 * \author Hannes Hauswedell <hannes.hauswedell AT fu-berlin.de>
 */

#pragma once

#include <seqan3/range/view/char_to.hpp>
#include <seqan3/range/view/concept.hpp>
#include <seqan3/range/view/complement.hpp>
#include <seqan3/range/view/convert.hpp>
#include <seqan3/range/view/rank_to.hpp>
#include <seqan3/range/view/to_char.hpp>
#include <seqan3/range/view/to_rank.hpp>
#include <seqan3/range/view/trim.hpp>

/*!\defgroup view View
 * \brief Views are "lazy range combinators" that offer modified views onto other ranges.
 * \ingroup range
 * \sa https://ericniebler.github.io/range-v3/index.html#range-views
 * \sa range/view/all.hpp
 *
 * SeqAn3 makes heavy use of views as defined in the
 * [Ranges Technical Specification](http://en.cppreference.com/w/cpp/experimental/ranges). Currently the
 * implementation is based on the [range-v3 library](https://github.com/ericniebler/range-v3) and all those views
 * are available in the namespace ranges::view, see
 * [the overview](https://ericniebler.github.io/range-v3/index.html#range-views) for more details.
 *
 * This submodule provides additional views, specifically for operations on biological data and
 * sequence analysis.
 *
 * \attention
 * To prevent naming conflicts, all SeqAn views are inside the namespace seqan3::view.
 *
 * \par Example
 *
 * ```cpp
 * dna4_vector vec{"ACGGTC"_dna4};
 * auto vec_view  = vec | view::complement;                         // == "TGCCAG" (but doesn't own any data)
 * auto vec_view2 = vec | ranges::view::reverse;                    // == "CTGGCA" (but doesn't own any data)
 *
 * // or in one line:
 * auto vec_view3 = vec | view::complement | ranges::view::reverse; // == "GACCGT" (but doesn't own any data)
 * ```
 */

/*!
 * \namespace seqan3::view
 * \brief The SeqAn3 namespace for views.
 *
 * Since views often have name clashes with regular functions and ranges they are implemented in the sub
 * namespace `view`.
 *
 * See the \link view view submodule \endlink of the range module for more details.
 */
