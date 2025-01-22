// SPDX-FileCopyrightText: 2006-2025 Knut Reinert & Freie Universität Berlin
// SPDX-FileCopyrightText: 2016-2025 Knut Reinert & MPI für molekulare Genetik
// SPDX-License-Identifier: BSD-3-Clause

/*!\file
 * \author Sara Hetzel <sara.hetzel AT fu-berlin.de>
 * \brief Genetic codes used for translating a triplet of nucleotides into an amino acid.
 */

#pragma once

#include <seqan3/core/platform.hpp>

namespace seqan3
{
/*!\brief Genetic codes used for translation of nucleotides into amino acids.
 *
 * \details
 * The numeric values of the enums correspond to the genbank transl_table values.
 * \sa https://www.ncbi.nlm.nih.gov/Taxonomy/Utils/wprintgc.cgi
 *
 * \experimentalapi{Experimental since version 3.1.}
 */
enum struct genetic_code : uint8_t
{
    canonical = 1, //!< [The Standard Code](https://www.ncbi.nlm.nih.gov/Taxonomy/Utils/wprintgc.cgi#SG1).
    //     vert_mitochondrial,
    //     yeast_mitochondrial,
    //     mold_mitochondrial,
    //     invert_mitochondrial,
    //     ciliate,
    //     flatworm_mitochondrial = 9,
    //     euplotid,
    //     prokaryote,
    //     alt_yeast,
    //     ascidian_mitochondrial,
    //     alt_flatworm_mitochondrial,
    //     blepharisma,
    //     chlorophycean_mitochondrial,
    //     trematode_mitochondrial = 21,
    //     scenedesmus_mitochondrial,
    //     thraustochytrium_mitochondrial,
    //     pterobranchia_mitochondrial,
    //     gracilibacteria
};
} // namespace seqan3
