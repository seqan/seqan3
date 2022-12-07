// -----------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2022, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2022, Knut Reinert & MPI für molekulare Genetik
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at: https://github.com/seqan/seqan3/blob/master/LICENSE.md
// -----------------------------------------------------------------------------------------------------

/*!\file
 * \author Christopher Pockrandt <christopher.pockrandt AT fu-berlin.de>
 * \brief Meta-header for the \link search Search module \endlink.
 */

/*!\defgroup search Search
 * \brief Data structures and approximate string search algorithms for large collection of text (e.g. DNA).
 *
 * \details
 *
 * Searching is a key component in many sequence analysis tools. The search module is a powerful and easy way to search
 * sequences in a large text or an arbitrary nested collection of texts. When it comes to searching, indices are a core
 * component for searching large amounts of data and are used for tools such as read mappers, assemblers or protein
 * search tools.
 *
 * \if KMER
 * SeqAn currently implements two kinds of indices: FM indices and k-mer indices (also known as q-gram indices).
 * Generally speaking, k-mer indices support very fast searching of exact k-mers (strings of length k) or k-mers with
 * predefined wildcard positions, which match with any character.
 * FM indices on the other hand are more versatile and work
 * with arbitrary pattern lengths and error numbers or positions.
 * \else
 * SeqAn currently implements only the FM index and a k-mer index is planned. The FM index works
 * with arbitrary pattern lengths and error numbers.
 * \endif
 *
 * \if DEV
 * \todo Elaborate on that (space consumption for growing k, maybe a rule of thumb).
 * \todo Add a description of the DREAM index here.
 * \endif
 *
 * # Search algorithm
 *
 * The Search module offers a unified search interface seqan3::search.
 * The function chooses the best search method based on the provided **index** and an optional **configuration**.
 *
 * \snippet test/snippet/search/search_with_user_callback.cpp search_with_config
 *
 * ### The search configuration
 *
 * The `seqan3::search` algorithm can be configured in multiple ways. You simply pass a respective `config` object to
 * the algorithm:
 *
 * The search configuration is independent of the index you use. It specifies for example how many errors are allowed
 * in the your search to still consider it successful.
 *
 * See detailed documentation on \ref search_configuration for details.
 *
 * ### The search result
 *
 * The result is a range over "hits". We call a *hit* the position of a successful search.
 * A *hit* in SeqAn is represented by the seqan3::search_result. You can configure what's part of the search result
 * with the \ref search_configuration_subsection_output.
 *
 * You can iterate over the results with a normal for loop:
 *
 * \snippet test/snippet/search/search.cpp search_result_loop
 *
 * ### Which index do I use?
 *
 * In order to choose the right Index for your use case, have a look at the next section
 * (\ref search_available_indices). The documentation of the respective index will lead you to a more detailed
 * description on how to use our search algorithm.
 *
 * # Available Indices \anchor search_available_indices
 *
 * ## seqan3::fm_index
 *
 * The seqan3::fm_index implements a original [FM Index](https://en.wikipedia.org/wiki/FM-index) with a trivial
 * backtracking approach. It is recommended for searches with no errors (exact searches). For eaxact searches, the
 * original FM index is slightly faster than the bidirectional FM Index and in general it is only half the size.
 *
 * ## seqan3::bi_fm_index
 *
 * The seqan3::bi_fm_index is a bidirectional FM index [1]. It improves the original FM index by allowing to extend the
 * query to the right and the left.  This make the bidirectional FM index more efficient than its predecessor when
 * searching with errors (approximate search).
 * The performance gain is enabled by using (optimum) *search schemes*.
 * Currently, using *search schemes* is only supported for searches with up to three errors.
 * In the future we plan to improve the search schemes to handle higher error counts.
 *
 * **Reference:**
 *
 * [1] Optimum Search Schemes for Approximate String Matching Using Bidirectional FM-Index. bioRxiv, 301085.
 * https://doi.org/10.1101/301085, *Kianfar, K., Pockrandt, C., Torkamandi, B., Luo, H., & Reinert, K. (2018).*
 *
 * \if DEV
 * ### Implementation details of Search Schemes
 *
 * A search scheme is a collection of searches, where each search `s` specifies the order in which the blocks are
 * searched (`s.pi`), the lower error bounds (`s.l`) and the upper error bounds (`s.u`). If the number of blocks that
 * the query sequences are split into are known at compile time, the data structure seqan3::detail::search is
 * recommended, otherwise one has to use seqan3::detail::search_dyn. The first one implements its member variables as an
 * a `std::array` of integers, the latter as a `std::vector` of integers.
 * Search search schemes are defined similiar. They are either implemented as a `std::array` of searches if the number
 * of searches is known at compile time, or as a `std::vector` if not.
 *
 * Precomputed optimum search schemes are represented as a `std::array` of seqan3::detail::search since both the
 * number of searches and the number of blocks are known at compile time. Search schemes computed at run time are
 * represented as `std::vector` of seqan3::detail::search_dyn.
 *
 * All optimum search schemes are disjoint, i.e. no error configuration is covered by more than one search. Thus there
 * is no need for a filtration phase after searching to remove duplicate hits when searching with hamming distance.
 * Allowing insertions and deletions will, however, lead to redundant reports of hits regardless of search schemes.
 *
 * \endif
 *
 */

#pragma once

#include <seqan3/search/configuration/all.hpp>
#include <seqan3/search/dream_index/all.hpp>
#include <seqan3/search/fm_index/all.hpp>
#include <seqan3/search/kmer_index/all.hpp>
#include <seqan3/search/search.hpp>
#include <seqan3/search/search_result.hpp>
#include <seqan3/search/views/all.hpp>
