// SPDX-FileCopyrightText: 2006-2025 Knut Reinert & Freie Universität Berlin
// SPDX-FileCopyrightText: 2016-2025 Knut Reinert & MPI für molekulare Genetik
// SPDX-License-Identifier: BSD-3-Clause

/*!\file
 * \author Hannes Hauswedell <hannes.hauswedell AT fu-berlin.de>
 * \brief Provides aliases for qualified.
 */

#pragma once

#include <seqan3/alphabet/concept.hpp>
#include <seqan3/alphabet/nucleotide/dna15.hpp>
#include <seqan3/alphabet/nucleotide/dna4.hpp>
#include <seqan3/alphabet/nucleotide/dna5.hpp>
#include <seqan3/alphabet/nucleotide/rna15.hpp>
#include <seqan3/alphabet/nucleotide/rna4.hpp>
#include <seqan3/alphabet/nucleotide/rna5.hpp>
#include <seqan3/alphabet/quality/concept.hpp>
#include <seqan3/alphabet/quality/phred42.hpp>
#include <seqan3/alphabet/quality/qualified.hpp>

namespace seqan3
{

/*!\brief An alphabet that stores a seqan3::dna4 letter and an seqan3::phred42 letter at each position.
 * \ingroup alphabet_quality
 * \details
 * \stableapi{Since version 3.1.}
 */
using dna4q = qualified<dna4, phred42>;

/*!\brief An alphabet that stores a seqan3::dna5 letter and an seqan3::phred42 letter at each position.
 * \ingroup alphabet_quality
 * \details
 * \stableapi{Since version 3.1.}
 */
using dna5q = qualified<dna5, phred42>;

/*!\brief An alphabet that stores a seqan3::rna4 letter and an seqan3::phred42 letter at each position.
 * \ingroup alphabet_quality
 * \details
 * \stableapi{Since version 3.1.}
 */
using rna4q = qualified<rna4, phred42>;

/*!\brief An alphabet that stores a seqan3::rna5 letter and an seqan3::phred42 letter at each position.
 * \ingroup alphabet_quality
 * \details
 * \stableapi{Since version 3.1.}
 */
using rna5q = qualified<rna5, phred42>;

/*!\brief An alphabet that stores a seqan3::dna15 letter and an seqan3::qualified letter at each position.
 * \ingroup alphabet_quality
 * \details
 * \stableapi{Since version 3.1.}
 */
using dna15q = qualified<dna15, phred42>;

/*!\brief An alphabet that stores a seqan3::rna15 letter and an seqan3::qualified letter at each position.
 * \ingroup alphabet_quality
 * \details
 * \stableapi{Since version 3.1.}
 */
using rna15q = qualified<rna15, phred42>;

} // namespace seqan3
