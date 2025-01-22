// SPDX-FileCopyrightText: 2006-2025 Knut Reinert & Freie Universität Berlin
// SPDX-FileCopyrightText: 2016-2025 Knut Reinert & MPI für molekulare Genetik
// SPDX-License-Identifier: BSD-3-Clause

/*!\file
 * \author Sara Hetzel <sara.hetzel AT fu-berlin.de>
 * \brief Provides functions for translating a triplet of nucleotides into an amino acid.
 */

#pragma once

#include <tuple>

#include <seqan3/alphabet/aminoacid/aa27.hpp>
#include <seqan3/alphabet/aminoacid/translation_details.hpp>
#include <seqan3/alphabet/aminoacid/translation_genetic_code.hpp>
#include <seqan3/core/range/type_traits.hpp>

namespace seqan3
{

// forwards:
class dna4;
class dna5;
class dna15;
class rna4;
class rna5;
class rna15;

/*!\brief Translate one nucleotide triplet into single amino acid (single nucleotide interface).
 * \ingroup alphabet_aminoacid
 * \tparam nucl_type The type of input nucleotides; must model seqan3::nucleotide_alphabet.
 * \param[in] n1 First nucleotide in triplet.
 * \param[in] n2 Second nucleotide in triplet.
 * \param[in] n3 Third nucleotide in triplet.
 *
 * \details
 *
 * Translates single nucleotides into amino acid according to given genetic code.
 *
 * ### Complexity
 *
 * Constant.
 *
 * ### Exceptions
 *
 * No-throw guarantee.
 *
 * \experimentalapi{Experimental since version 3.1.}
 */
template <genetic_code gc = genetic_code::canonical, nucleotide_alphabet nucl_type>
constexpr aa27 translate_triplet(nucl_type const & n1, nucl_type const & n2, nucl_type const & n3) noexcept
{
    if constexpr (std::same_as<nucl_type, dna4> || std::same_as<nucl_type, dna5> || std::same_as<nucl_type, dna15>)
    {
        // table exists for dna15 and is generated for dna4 and dna5 (compile time ok, because small)
        return seqan3::detail::translation_table<nucl_type, gc>::value[to_rank(n1)][to_rank(n2)][to_rank(n3)];
    }
    else if constexpr (std::same_as<nucl_type, rna4> || std::same_as<nucl_type, rna5> || std::same_as<nucl_type, rna15>)
    {
        using rna2dna_t =
            std::conditional_t<std::same_as<nucl_type, rna4>,
                               dna4,
                               std::conditional_t<std::same_as<nucl_type, rna5>,
                                                  dna5,
                                                  std::conditional_t<std::same_as<nucl_type, rna15>, dna15, void>>>;

        // we can use dna's tables, because ranks are identical
        return seqan3::detail::translation_table<rna2dna_t, gc>::value[to_rank(n1)][to_rank(n2)][to_rank(n3)];
    }
    else // composites or user defined nucleotide
    {
        // we cast to dna15; slightly slower run-time, but lot's of compile time saved for large alphabets.
        // (nucleotide types can be converted to dna15 by definition)
        return seqan3::detail::translation_table<dna15, gc>::value[to_rank(static_cast<dna15>(n1))][to_rank(
            static_cast<dna15>(n2))][to_rank(static_cast<dna15>(n3))];
    }
}

} // namespace seqan3
