// -----------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2019, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2019, Knut Reinert & MPI für molekulare Genetik
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at: https://github.com/seqan/seqan3/blob/master/LICENSE.md
// -----------------------------------------------------------------------------------------------------

 /*!\file
  * \author Christopher Pockrandt <christopher.pockrandt AT fu-berlin.de>
  * \brief Meta-header for the Search Algorithm module.
  *
  * \defgroup submodule_search_algorithm Algorithm
  * \brief Provides seqan3::search.
  * \ingroup search
  *
  * \details
  * # Search Algorithms
  *
  * The Search module offers a simple unified interface that allows searching FM indices and choosing the best
  * algorithm based on the index at hand.
  *
  * # FM Indices
  *
  * The search algorithms for FM indices implement either a trivial backtracking approach or an optimum search scheme.
  * Latter are currently only available for searches with up to and including three errors using bidirectional indices.
  * The optimum search schemes will be improved in the future to handle unidirectional indices and higher error counts.
  *
  * \if DEV
  *
  * ## Implementation details of Search Schemes
  *
  * A search scheme is a collection of searches, where each search `s` specifies the order in which the blocks are
  * searched (`s.pi`), the lower error bounds (`s.l`) and the upper error bounds (`s.u`). If the number of blocks that
  * the query sequences are split into are known at compile time, the data structure `search` is recommended, otherwise
  * one has to use seqan3::detail::search_dyn. The first one implements its member variables as an an `std::array` of
  * integers, the latter as an `std::vector` of integers.
  * Similarily search schemes are defined. They are either implemented as an `std::array` of searches if the number of
  * searches is known at compile time, or as an `std::vector` if not.
  *
  * Precomputed optimum search schemes are represented as an `std::array` of seqan3::detail::search since both the
  * number of searches and the number of blocks are known at compile time. Search schemes computed at run time are
  * represented as `std::vector` of seqan3::detail::search_dyn.
  *
  * All optimum search schemes are disjoint, i.e. no error configuration is covered by more than one search. Thus there
  * is no need for a filtration phase after searching to remove duplicate hits. Allowing insertions and deletions will
  * lead to redundant reports of hits regardless of search schemes.
  *
  * \endif
  *
  * ## Reference
  *
  * Kianfar, K., Pockrandt, C., Torkamandi, B., Luo, H., & Reinert, K. (2018).
  *
  * Optimum Search Schemes for Approximate String Matching Using Bidirectional FM-Index. bioRxiv, 301085.
  * https://doi.org/10.1101/301085
  *
  * \if DEV
  * # K-mer Indices
  *
  * \todo Rewrite landing page.
  * \endif
  *
  */

#pragma once

#include <seqan3/search/algorithm/configuration/all.hpp>
#include <seqan3/search/algorithm/search.hpp>
