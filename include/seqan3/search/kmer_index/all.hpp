// -----------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2020, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2020, Knut Reinert & MPI für molekulare Genetik
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at: https://github.com/seqan/seqan3/blob/master/LICENSE.md
// -----------------------------------------------------------------------------------------------------

/*!\file
 * \author Enrico Seiler <enrico.seiler AT fu-berlin.de>
 * \brief Meta-header for the kmer_index module.
 *
 * \defgroup submodule_kmer_index k-mer Index
 * \ingroup search
 * \brief Implementation of shapes for a k-mer Index.
 *
 * \note The k-mer index is not yet implemented.
 *
 * \details
 *
 * A k-mer index is a data structure that stores all occurrences of k-mers in a text.
 * A k-mer can be either an exact string of length k or it can contain one or more wildcards,
 * which denote positions of arbitrary characters.
 *
 * The index is very fast for retrieving all occurrences of a k-mer in the underlying text.
 * Usually the query length (k) is small and the underlying text is very large.
 * The parameter k and the position(s) of wildcards must be fixed at index creation with
 * seqan3::ungapped or seqan3::shape.
 */

#pragma once

#include <seqan3/search/kmer_index/shape.hpp>
