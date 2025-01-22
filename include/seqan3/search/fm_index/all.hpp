// SPDX-FileCopyrightText: 2006-2025 Knut Reinert & Freie Universität Berlin
// SPDX-FileCopyrightText: 2016-2025 Knut Reinert & MPI für molekulare Genetik
// SPDX-License-Identifier: BSD-3-Clause

/*!\file
 * \author Christopher Pockrandt <christopher.pockrandt AT fu-berlin.de>
 * \brief Meta-header for the \link search_fm_index Search / FM Index submodule \endlink.
 */

/*!\defgroup search_fm_index FM Index
 * \brief Provides seqan3::fm_index and seqan3::bi_fm_index as well as respective cursors.
 * \ingroup search
 *
 * \details
 * # FM Indices
 *
 * FM indices are text indices similar to suffix trees or suffix arrays which are based on the Burrow Wheeler
 * transform and a sampled suffix array. FM indices are significantly smaller in space without sacrificing speed. They
 * also allow for adjusting the speed respectively space by choosing the underlying data structures accordingly or
 * modyfing the sampling rate of the sampled suffix array.
 *
 * The FM indices are based on the <a href="https://github.com/xxsds/sdsl-lite">SDSL 3</a> (succinct data structure
 * library). You are able to specify the underlying implementation of the SDSL to adjust it to your needs as well as
 * choose one of the preconfigured indices that are suitable for common applications in sequence analysis.
 *
 * For technical reasons you can currently only build indices over a seqan3::alphabet if its
 * seqan3::alphabet_size is smaller or equal 256.
 *
 * You can choose between unidirectional and bidirectional FM indices (which can be thought of suffix trees
 * and affix trees, i.e. a combination of suffix and prefix trees being able to search a pattern from left to
 * right, right to left and character by character in any arbitrary order). Roughly speaking bidirectional
 * FM indices are more powerful for approximate string matching at the cost of a higher space consumption
 * (between a factor of 5 and 9 of the input size depending on the configuration).
 *
 * # FM Index Cursors
 *
 * Index Cursors are lightweight objects, i.e. they are cheap to copy.
 *
 * Note that the SeqAn index cursor, although having similar behaviour, doesn't model any of the standard library
 * iterator concepts, not even std::input_or_output_iterator.
 */

#pragma once

#include <seqan3/search/fm_index/bi_fm_index.hpp>
#include <seqan3/search/fm_index/bi_fm_index_cursor.hpp>
#include <seqan3/search/fm_index/concept.hpp>
#include <seqan3/search/fm_index/fm_index.hpp>
#include <seqan3/search/fm_index/fm_index_cursor.hpp>
