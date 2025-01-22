// SPDX-FileCopyrightText: 2006-2025 Knut Reinert & Freie Universität Berlin
// SPDX-FileCopyrightText: 2016-2025 Knut Reinert & MPI für molekulare Genetik
// SPDX-License-Identifier: BSD-3-Clause

/*!\file
 * \author Joerg Winkler <j.winkler AT fu-berlin.de>
 * \brief Meta-header for the \link alphabet_structure Alphabet / Structure submodule \endlink.
 */

/*!\defgroup alphabet_structure Structure
 * \brief Provides types to represent single elements of RNA and protein structures.
 * \ingroup alphabet
 * \see alphabet
 *
 * \details
 *
 * The following alphabets are currently supported in SeqAn. Please see the format's page for more details.
 *
 * Name                                     | Chars               | Description
 * ---------------------------------------- | ------------------------ | -----------
 * [Dot Bracket](@ref seqan3::dot_bracket3) | `().`                    | Simple annotation that defines base pairs. No pseudoknots allowed.
 * [WUSS](@ref seqan3::wuss)                | `.<>:,-_~;()[]{}AaBb...` | Annotation that provides further markups and pseudoknots.
 * [DSSP](@ref seqan3::dssp9)               | `HBEGITSCX`              | Structure encoding for proteins.
 */

#pragma once

#include <seqan3/alphabet/structure/concept.hpp>
#include <seqan3/alphabet/structure/dot_bracket3.hpp>
#include <seqan3/alphabet/structure/dssp9.hpp>
#include <seqan3/alphabet/structure/structured_aa.hpp>
#include <seqan3/alphabet/structure/structured_rna.hpp>
#include <seqan3/alphabet/structure/wuss.hpp>
