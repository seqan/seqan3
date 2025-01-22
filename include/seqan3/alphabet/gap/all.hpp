// SPDX-FileCopyrightText: 2006-2025 Knut Reinert & Freie Universität Berlin
// SPDX-FileCopyrightText: 2016-2025 Knut Reinert & MPI für molekulare Genetik
// SPDX-License-Identifier: BSD-3-Clause

/*!\file
 * \author Marcel Ehrhardt <marcel.ehrhardt AT fu-berlin.de>
 * \brief Meta-header for the \link alphabet_gap Alphabet / Gap submodule \endlink.
 */

/*!\defgroup alphabet_gap Gap
 * \brief Provides the gap alphabet and functionality to make an alphabet a gapped alphabet.
 * \ingroup alphabet
 * \see alphabet
 *
 * \details
 *
 * ### Introduction
 *
 * The gap symbol (`-`) is used in alignments to represent an interruption in an alignment, usually the result of an
 * insertion or deletion. The seqan3::gap alphabet represents this (single) gap symbol and satisfies the
 * seqan3::alphabet.
 *
 * The main purpose of seqan3::gap is to be combined with other alphabets. This can easily be achieved by using the
 * seqan3::gapped<> template which transforms any other alphabet to be a composite of that alphabet + the gap
 * character.
 */

#pragma once

#include <seqan3/alphabet/gap/gap.hpp>
#include <seqan3/alphabet/gap/gapped.hpp>
