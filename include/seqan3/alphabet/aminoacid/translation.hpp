// -----------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2020, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2020, Knut Reinert & MPI für molekulare Genetik
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at: https://github.com/seqan/seqan3/blob/master/LICENSE.md
// -----------------------------------------------------------------------------------------------------

/*!\file
 * \author Sara Hetzel <sara.hetzel AT fu-berlin.de>
 * \brief Provides functions for translating a triplet of nucleotides into an amino acid.
 */

#pragma once

#include <tuple>

#include <seqan3/alphabet/aminoacid/aa27.hpp>
#include <seqan3/alphabet/aminoacid/translation_genetic_code.hpp>
#include <seqan3/alphabet/aminoacid/translation_details.hpp>
#include <seqan3/core/type_traits/range.hpp>

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
 * \ingroup aminoacid
 * \tparam nucl_type The type of input nucleotides.
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
 */
template <genetic_code gc = genetic_code::CANONICAL, nucleotide_alphabet nucl_type>
constexpr aa27 translate_triplet(nucl_type const & n1, nucl_type const & n2, nucl_type const & n3) noexcept
{
    if constexpr (std::same_as<nucl_type, dna4> || std::same_as<nucl_type, dna5> || std::same_as<nucl_type, dna15>)
    {
        // table exists for dna15 and is generated for dna4 and dna5 (compile time ok, because small)
        return seqan3::detail::translation_table<nucl_type, gc>::VALUE[to_rank(n1)][to_rank(n2)][to_rank(n3)];
    }
    else if constexpr (std::same_as<nucl_type, rna4> || std::same_as<nucl_type, rna5> || std::same_as<nucl_type, rna15>)
    {
        using rna2dna_t = std::conditional_t<std::same_as<nucl_type, rna4>,  dna4,
                          std::conditional_t<std::same_as<nucl_type, rna5>,  dna5,
                          std::conditional_t<std::same_as<nucl_type, rna15>, dna15, void>>>;

        // we can use dna's tables, because ranks are identical
        return seqan3::detail::translation_table<rna2dna_t, gc>::VALUE[to_rank(n1)][to_rank(n2)][to_rank(n3)];
    }
    else // composites or user defined nucleotide
    {
        // we cast to dna15; slightly slower run-time, but lot's of compile time saved for large alphabets.
        // (nucleotide types can be converted to dna15 by definition)
        return seqan3::detail::translation_table<dna15, gc>::VALUE[to_rank(static_cast<dna15>(n1))]
                                                                  [to_rank(static_cast<dna15>(n2))]
                                                                  [to_rank(static_cast<dna15>(n3))];
    }
}

/*!\brief Translate one nucleotide triplet into single amino acid (tuple interface).
 * \ingroup aminoacid
 * \tparam tuple_type Type of `input_tuple`. Usually std::tuple, but similar types like std::array
 * and seqan3::pod_tuple are also supported.
 * \param[in] input_tuple Triplet of nucleotides that should be converted to amino acid.
 * \details
 *
 * Translates std::tuple or std::array with 3 nucleotides into amino acid according to given genetic code.
 *
 * ### Complexity
 *
 * Constant.
 *
 * ### Exceptions
 *
 * No-throw guarantee.
 *
 * \deprecated Use seqan3::translate_triplet(nucl_type const & n1, nucl_type const & n2, nucl_type const & n3) instead.
 */
template <genetic_code gc = genetic_code::CANONICAL, typename tuple_type>
//!\cond
    requires (std::tuple_size<tuple_type>::value == 3) &&
             nucleotide_alphabet<std::tuple_element_t<0, tuple_type>> &&
             nucleotide_alphabet<std::tuple_element_t<1, tuple_type>> &&
             nucleotide_alphabet<std::tuple_element_t<2, tuple_type>>
//!\endcond
constexpr aa27 translate_triplet SEQAN3_DEPRECATED_310 (tuple_type const & input_tuple) noexcept
{
    return translate_triplet(std::get<0>(input_tuple), std::get<1>(input_tuple), std::get<2>(input_tuple));
}

/*!\brief Translate one nucleotide triplet into single amino acid (range interface).
 * \ingroup aminoacid
 * \tparam range_type Type of input_range; must satisfy std::ranges::input_range.
 * \param[in] input_range Range of three nucleotides that should be converted to amino acid.
 *
 * \details
 *
 * Translates range with 3 nucleotides into amino acid according to given genetic code.
 *
 * ### Complexity
 *
 * Constant.
 *
 * ### Exceptions
 *
 * Strong exception guarantee (never modifies data).
 *
 * \deprecated Use seqan3::translate_triplet(nucl_type const & n1, nucl_type const & n2, nucl_type const & n3) instead.
 */
template <genetic_code gc = genetic_code::CANONICAL, std::ranges::input_range range_type>
//!\cond
    requires nucleotide_alphabet<std::ranges::range_reference_t<std::decay_t<range_type>>>
//!\endcond
constexpr aa27 translate_triplet SEQAN3_DEPRECATED_310 (range_type && input_range)
{
    auto n1 = std::ranges::begin(input_range);
    auto n2 = ++n1;
    auto n3 = ++n2;

    assert(n1 != std::ranges::end(input_range));
    assert(n2 != std::ranges::end(input_range));
    assert(n3 != std::ranges::end(input_range));

    return translate_triplet(*n1, *n2, *n3);
}

/*!\brief Translate one nucleotide triplet into single amino acid (range interface, input range allows random access).
 * \ingroup aminoacid
 * \tparam rng_t Type of input_range; must satisfy std::ranges::random_access_range.
 * \param[in] input_range Range of three nucleotides that should be converted to amino acid.
 *
 * \details
 *
 * Translates range with 3 nucleotides into amino acid according to given genetic code.
 *
 * ### Complexity
 *
 * Constant.
 *
 * ### Exceptions
 *
 * Strong exception guarantee (never modifies data).
 *
 * \deprecated Use seqan3::translate_triplet(nucl_type const & n1, nucl_type const & n2, nucl_type const & n3) instead.
 */
template <genetic_code gc = genetic_code::CANONICAL, std::ranges::random_access_range rng_t>
//!\cond
    requires nucleotide_alphabet<std::ranges::range_reference_t<std::decay_t<rng_t>>>
//!\endcond
constexpr aa27 translate_triplet SEQAN3_DEPRECATED_310 (rng_t && input_range)
{
    assert(std::ranges::begin(input_range) != std::ranges::end(input_range));
    assert(std::ranges::begin(input_range) + 1 != std::ranges::end(input_range));
    assert(std::ranges::begin(input_range) + 2 != std::ranges::end(input_range));

    return translate_triplet(input_range[0], input_range[1], input_range[2]);
}

} // namespace seqan3
