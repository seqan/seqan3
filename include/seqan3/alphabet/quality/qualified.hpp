// -----------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2019, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2019, Knut Reinert & MPI für molekulare Genetik
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at: https://github.com/seqan/seqan3/blob/master/LICENSE
// -----------------------------------------------------------------------------------------------------

/*!\file
 * \author Hannes Hauswedell <hannes.hauswedell AT fu-berlin.de>
 * \brief Contains quality alphabet compositions.
 */

#pragma once

#include <iostream>
#include <string>
#include <utility>

#include <seqan3/alphabet/composition/cartesian_composition.hpp>
#include <seqan3/alphabet/nucleotide/concept.hpp>

namespace seqan3
{

/*!\brief Joins an arbitrary alphabet with a quality alphabet.
 * \ingroup quality
 * \tparam sequence_alphabet_t Type of the alphabet; must satisfy seqan3::alphabet_concept.
 * \tparam quality_alphabet_t  Type of the quality; must satisfy seqan3::quality_concept.
 * \implements seqan3::quality_concept
 * \implements seqan3::detail::constexpr_alphabet_concept
 * \implements seqan3::trivially_copyable_concept
 * \implements seqan3::standard_layout_concept
 *
 *
 * This composition pairs an arbitrary alphabet with a quality alphabet, where
 * each alphabet character is stored together with its quality score in a
 * single value. That way, you can can conveniently access the character and
 * score information at each position of the qualified-sequence.
 * The use case that this was designed for is a nucleotide sequence with
 * corresponding quality scores, e.g. obtained when reading in a FASTQ file
 * of Illumina reads.
 * The composition also allows to store quality scores for different or extended
 * alphabets like a `qualified<char, phred42>` or a `qualified<gapped<dna4>, phred42>`
 * sequence.
 * The rank values correspond to numeric values in the size of the composition,
 * while the character values are taken from the sequence alphabet and the phred
 * values are taken from the quality alphabet.
 *
 * As with all `seqan3::cartesian_composition` s you may access the individual
 * alphabet letters in regular c++ tuple notation, i.e. `get<0>(t)` and objects
 * can be brace-initialised with the individual members.
 *
 * \snippet test/snippet/alphabet/quality/qualified.cpp general
 *
 * This seqan3::cartesian_composition itself fulfils both seqan3::alphabet_concept and seqan3::quality_concept.
 */
template <alphabet_concept sequence_alphabet_t, quality_concept quality_alphabet_t>
class qualified :
    public cartesian_composition<qualified<sequence_alphabet_t, quality_alphabet_t>,
                                 sequence_alphabet_t, quality_alphabet_t>
{
private:
    //!\brief The base type.
    using base_type = cartesian_composition<qualified<sequence_alphabet_t, quality_alphabet_t>,
                                            sequence_alphabet_t, quality_alphabet_t>;

public:
    //!\brief First template parameter as member type.
    using sequence_alphabet_type = sequence_alphabet_t;
    //!\brief Second template parameter as member type.
    using quality_alphabet_type = quality_alphabet_t;

    //!\brief Equals the char_type of sequence_alphabet_type.
    using char_type = underlying_char_t<sequence_alphabet_type>;
    //!\brief Equals the phred_type of the quality_alphabet_type.
    using phred_type = underlying_phred_t<quality_alphabet_type>;

    /*!\name Constructors, destructor and assignment
     * \{
     */
    constexpr qualified() noexcept : base_type{} {}
    constexpr qualified(qualified const &) = default;
    constexpr qualified(qualified &&) = default;
    constexpr qualified & operator =(qualified const &) = default;
    constexpr qualified & operator =(qualified &&) = default;
    ~qualified() = default;

    using base_type::base_type; // Inherit non-default constructors

    // Inherit operators from base
    using base_type::operator=;
    using base_type::operator==;
    using base_type::operator!=;
    using base_type::operator>=;
    using base_type::operator<=;
    using base_type::operator<;
    using base_type::operator>;

    //!\copydoc cartesian_composition::cartesian_composition(component_type const alph)
    SEQAN3_DOXYGEN_ONLY(( constexpr qualified(component_type const alph) noexcept {} ))
    //!\copydoc cartesian_composition::cartesian_composition(indirect_component_type const alph)
    SEQAN3_DOXYGEN_ONLY(( constexpr qualified(indirect_component_type const alph) noexcept {} ))
    //!\copydoc cartesian_composition::operator=(component_type const alph)
    SEQAN3_DOXYGEN_ONLY(( constexpr qualified & operator=(component_type const alph) noexcept {} ))
    //!\copydoc cartesian_composition::operator=(indirect_component_type const alph)
    SEQAN3_DOXYGEN_ONLY(( constexpr qualified & operator=(indirect_component_type const alph) noexcept {} ))
    //!\}

    /*!\name Write functions
     * \{
     */
    //!\brief Assign from a character. This modifies the internal sequence letter.
    constexpr qualified & assign_char(char_type const c) noexcept
    {
        seqan3::assign_char(get<0>(*this), c);
        return *this;
    }

    //!\brief Strict assign from a character. This modifies the internal sequence letter.
    qualified & assign_char_strict(char_type const c)
    {
        seqan3::assign_char_strict(get<0>(*this), c);
        return *this;
    }

    //!\brief Assign from a phred value. This modifies the internal quality letter.
    constexpr qualified & assign_phred(phred_type const c) noexcept
    {
        seqan3::assign_phred(get<1>(*this), c);
        return *this;
    }
    //!\}

    /*!\name Read functions
     * \{
     */
    //!\brief Return the phred value. This reads the internal quality letter.
    constexpr phred_type to_phred() const noexcept
    {
        return seqan3::to_phred(get<1>(*this));
    }

    //!\brief Return a character. This reads the internal sequence letter.
    constexpr char_type to_char() const noexcept
    {
        return seqan3::to_char(get<0>(*this));
    }

    /*!\brief Return a qualified where the quality is preserved, but the sequence letter is complemented.
     * \sa seqan3::complement
     * \sa seqan3::nucleotide_concept::complement
     */
    constexpr qualified complement() const noexcept
        requires nucleotide_concept<sequence_alphabet_t>
    {
        using seqan3::complement;
        return qualified{complement(get<0>(*this)), get<1>(*this)};
    }
    //!\}

    //!\brief Validate whether a character is valid in the sequence alphabet.
    static constexpr bool char_is_valid(char_type const c) noexcept
    {
        return char_is_valid_for<sequence_alphabet_type>(c);
    }
};

//!\brief Type deduction guide enables usage of qualified without specifying template args.
//!\relates qualified
template <typename sequence_alphabet_type, typename quality_alphabet_type>
qualified(sequence_alphabet_type &&, quality_alphabet_type &&)
    -> qualified<std::decay_t<sequence_alphabet_type>, std::decay_t<quality_alphabet_type>>;

} // namespace seqan3
