// -----------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2019, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2019, Knut Reinert & MPI für molekulare Genetik
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at: https://github.com/seqan/seqan3/blob/master/LICENSE.md
// -----------------------------------------------------------------------------------------------------

/*!\file
 * \author Christopher Pockrandt <christopher.pockrandt AT fu-berlin.de>
 * \brief Meta-header for the search module.
 *
 * \defgroup search Search
 * \brief Contains datastructures and algorithms for searching.
 *
 * \details
 * # Introduction
 *
 * Searching is a key component in many sequence analysis tools. The search module is a powerful and easy way to search
 * sequences in a large text or an arbitrary nested collection of texts. When it comes to searching, indices are a core
 * component for searching large amounts of data and are used for tools such as read mappers, assemblers or protein
 * search tools.
 *
 * \if DEV
 * There are currently two major kind of indices: FM indices and k-mer indices (also known as q-gram
 * indices).
 *
 * \todo Elaborate on that (space consumption for growing k, maybe a rule of thumb).
 *
 * Generally speaking k-mer indices support very fast searching of exact k-mers (strings of length k) or k-mers with
 * predefined wildcard positions that do not have to match. FM indices on the other hand are more versatile and work
 * with arbitrary pattern lengths and error numbers / positions.
 * \todo k-mer indices are coming soon. Stay tuned!
 * \endif
 *
 * SeqAn3 currently supports very fast FM indices. For more details visit the \ref submodule_fm_index
 * "FM Index submodule".
 *
 */

#pragma once

#include <seqan3/search/algorithm/all.hpp>
#include <seqan3/search/fm_index/all.hpp>
