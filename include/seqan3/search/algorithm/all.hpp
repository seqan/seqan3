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
  * \author Christopher Pockrandt <christopher.pockrandt AT fu-berlin.de>
  * \brief Meta-header for the Search Algorithm module.
  *
  * \defgroup submodule_search_algorithm Algorithm
  * \ingroup search
  *
  * \details
  * ## Search Algorithms
  *
  * The Search module offers a simple unified interface that allows searching SeqAn3 indices such as FM indices or k-mer
  * indices and choosing the best algorithm based on the index at hand.
  *
  * ## FM Indices
  *
  * The search algorithms for FM indices implement either a trivial backtracking approach or an optimum search scheme.
  * Latter are currently only available for searches with up to and including three errors using bidirectional indices.
  * The optimum search schemes will be improved in the future to handle unidirectional indices and higher error counts.
  *
  * ### Implementation details of Search Schemes
  *
  * \todo Remove this from user documentation (cond / endcond does not work within doxygen comment)
  *
  * A search scheme is a collection of searches, where each search `s` specifies the order in which the blocks are
  * searched (`s.pi`), the lower error bounds (`s.l`) and the upper error bounds (`s.u`). If the number of blocks that
  * the query sequences are split into are known at compile time, the data structure `search` is recommended, otherwise
  * one has to use `search_dyn`. The first one implements its member variables as an an `std::array` of integers, the
  * latter as an `std::vector` of integers.
  * Similarily search schemes are defined. They are either implemented as an `std::array` of searches if the number of
  * searches is known at compile time, or as an `std::vector` if not.
  *
  * Precomputed optimum search schemes are represented as an `std::array` of `search` since both the number of searches
  * and the number of blocks are known at compile time. Search schemes computed at running time are represented as
  * `std::vector` of `search_dyn`.
  *
  * All optimum search schemes are disjoint, i.e. no error configuration is covered by more than one search. Thus there
  * is no need for a filteration phase after searching to remove duplicate hits. This applies only when allowing for
  * substitutions only as for insertions and deletions lead to reundant reportings of hits regardless of search schemes.
  *
  * ### Reference
  *
  * Kianfar, K., Pockrandt, C., Torkamandi, B., Luo, H., & Reinert, K. (2018).
  *
  * Optimum Search Schemes for Approximate String Matching Using Bidirectional FM-Index. bioRxiv, 301085.
  *
  * ## K-mer Indices
  *
  * \todo Rewrite landing page.
  *
  */

#pragma once

#include <seqan3/search/algorithm/configuration/all.hpp>
#include <seqan3/search/algorithm/search.hpp>
