// -----------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2020, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2020, Knut Reinert & MPI für molekulare Genetik
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at: https://github.com/seqan/seqan3/blob/master/LICENSE.md
// -----------------------------------------------------------------------------------------------------

/*!\file
 * \author Christopher Pockrandt <christopher.pockrandt AT fu-berlin.de>
 * \brief Meta-header for the search module.
 *
 * \defgroup search Search
 * \brief Data structures and algorithms for the search of query sequences in a large collection of text.
 *
 * \details
 * # Introduction
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
 * The Search module offers a simple unified interface for searching a query in a large indexed text.
 * The algorithm chooses the best search method based on the provided index.
 *
 * # FM index
 *
 * The search algorithms for FM indices implement either a trivial backtracking approach or an optimum search scheme.
 * The latter are currently only available for searches with up to three errors using bidirectional indices.
 * In the future we plan to improve the optimum search schemes to handle higher error counts.
 *
 * \if DEV
 * ### Implementation details of Search Schemes
 *
 * A search scheme is a collection of searches, where each search `s` specifies the order in which the blocks are
 * searched (`s.pi`), the lower error bounds (`s.l`) and the upper error bounds (`s.u`). If the number of blocks that
 * the query sequences are split into are known at compile time, the data structure seqan3::detail::search is
 * recommended, otherwise one has to use seqan3::detail::search_dyn. The first one implements its member variables as an
 * an `std::array` of integers, the latter as an `std::vector` of integers.
 * Search search schemes are defined similiar. They are either implemented as an `std::array` of searches if the number
 * of searches is known at compile time, or as an `std::vector` if not.
 *
 * Precomputed optimum search schemes are represented as an `std::array` of seqan3::detail::search since both the
 * number of searches and the number of blocks are known at compile time. Search schemes computed at run time are
 * represented as `std::vector` of seqan3::detail::search_dyn.
 *
 * All optimum search schemes are disjoint, i.e. no error configuration is covered by more than one search. Thus there
 * is no need for a filtration phase after searching to remove duplicate hits when searching with hamming distance.
 * Allowing insertions and deletions will, however, lead to redundant reports of hits regardless of search schemes.
 *
 * \endif
 *
 * **Reference:**
 *
 * Kianfar, K., Pockrandt, C., Torkamandi, B., Luo, H., & Reinert, K. (2018).
 *
 * Optimum Search Schemes for Approximate String Matching Using Bidirectional FM-Index. bioRxiv, 301085.
 * https://doi.org/10.1101/301085
 *
 *
 * # K-mer index
 *
 * A k-mer index can be used to efficiently retrieve all occurrences of a certain k-mer in the text.
 * The k-mer can be either an exact string of length k or it can contain one or more wildcards,
 * which denote positions of arbitrary characters.
 *
 * An exact k-mer is represented as seqan3::ungapped, and wildcards can be defined with seqan3::shape.
 * Please check the respective documentation for details and examples.
 *
 * \note The k-mer index is not yet implemented.
 * \sa seqan3::views::kmer_hash
 * \sa seqan3::views::minimiser
 */

#pragma once

#include <seqan3/search/configuration/all.hpp>
#include <seqan3/search/fm_index/all.hpp>
#include <seqan3/search/search.hpp>
