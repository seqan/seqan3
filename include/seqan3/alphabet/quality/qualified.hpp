// SPDX-FileCopyrightText: 2006-2025 Knut Reinert & Freie Universität Berlin
// SPDX-FileCopyrightText: 2016-2025 Knut Reinert & MPI für molekulare Genetik
// SPDX-License-Identifier: BSD-3-Clause

/*!\file
 * \author Hannes Hauswedell <hannes.hauswedell AT fu-berlin.de>
 * \brief Provides quality alphabet composites.
 */

#pragma once

#include <seqan3/alphabet/composite/alphabet_tuple_base.hpp>
#include <seqan3/alphabet/nucleotide/concept.hpp>
#include <seqan3/alphabet/quality/concept.hpp>

namespace seqan3
{

/*!\brief Joins an arbitrary alphabet with a quality alphabet.
 * \ingroup alphabet_quality
 * \tparam sequence_alphabet_t Type of the alphabet; must satisfy seqan3::writable_alphabet.
 * \tparam quality_alphabet_t  Type of the quality; must satisfy seqan3::writable_quality_alphabet.
 * \implements seqan3::writable_quality_alphabet
 * \implements seqan3::detail::writable_constexpr_alphabet
 * \implements seqan3::trivially_copyable
 * \implements seqan3::standard_layout
 * \implements std::regular
 *
 *
 * This composite pairs an arbitrary alphabet with a quality alphabet, where
 * each alphabet character is stored together with its quality score in a
 * single value. That way, you can can conveniently access the character and
 * score information at each position of the qualified-sequence.
 * The use case that this was designed for is a nucleotide sequence with
 * corresponding quality scores, e.g. obtained when reading in a FASTQ file
 * of Illumina reads.
 * The composite also allows to store quality scores for different or extended
 * alphabets like a `qualified<char, phred42>` or a `qualified<gapped<dna4>, phred42>`
 * sequence.
 * The rank values correspond to numeric values in the size of the composite,
 * while the character values are taken from the sequence alphabet and the Phred score
 * values are taken from the quality alphabet.
 *
 * As with all `seqan3::alphabet_tuple_base` s you may access the individual
 * alphabet letters in regular c++ tuple notation, i.e. `get<0>(t)` and objects
 * can be brace-initialised with the individual members.
 *
 * \include test/snippet/alphabet/quality/qualified.cpp
 *
 * This seqan3::alphabet_tuple_base itself fulfils both seqan3::writable_alphabet and seqan3::writable_quality_alphabet.
 *
 * \stableapi{Since version 3.1.}
 */
template <writable_alphabet sequence_alphabet_t, writable_quality_alphabet quality_alphabet_t>
class qualified :
    public alphabet_tuple_base<qualified<sequence_alphabet_t, quality_alphabet_t>,
                               sequence_alphabet_t,
                               quality_alphabet_t>
{
private:
    //!\brief The base type.
    using base_type = alphabet_tuple_base<qualified<sequence_alphabet_t, quality_alphabet_t>,
                                          sequence_alphabet_t,
                                          quality_alphabet_t>;

public:
    /*!\brief First template parameter as member type.
     * \details
     * \stableapi{Since version 3.1.}
     */
    using sequence_alphabet_type = sequence_alphabet_t;
    /*!\brief Second template parameter as member type.
     * \details
     * \stableapi{Since version 3.1.}
     */
    using quality_alphabet_type = quality_alphabet_t;

    /*!\brief Equals the char_type of sequence_alphabet_type.
     * \details
     * \stableapi{Since version 3.1.}
     */
    using char_type = alphabet_char_t<sequence_alphabet_type>;
    /*!\brief Equals the phred_type of the quality_alphabet_type.
     * \details
     * \stableapi{Since version 3.1.}
     */
    using phred_type = alphabet_phred_t<quality_alphabet_type>;

    /*!\name Constructors, destructor and assignment
     * \{
     */
    constexpr qualified() noexcept = default;                              //!< Defaulted.
    constexpr qualified(qualified const &) noexcept = default;             //!< Defaulted.
    constexpr qualified(qualified &&) noexcept = default;                  //!< Defaulted.
    constexpr qualified & operator=(qualified const &) noexcept = default; //!< Defaulted.
    constexpr qualified & operator=(qualified &&) noexcept = default;      //!< Defaulted.
    ~qualified() noexcept = default;                                       //!< Defaulted.

    // Inherit from base:
    using base_type::alphabet_size;
    using base_type::base_type; // non-default constructors
    using base_type::to_rank;
    using base_type::operator=;

    //!\copydoc alphabet_tuple_base::alphabet_tuple_base(component_type const alph)
    SEQAN3_DOXYGEN_ONLY((constexpr qualified(component_type const alph) noexcept {}))
    //!\copydoc alphabet_tuple_base::alphabet_tuple_base(indirect_component_type const alph)
    SEQAN3_DOXYGEN_ONLY((constexpr qualified(indirect_component_type const alph) noexcept {}))
    //!\copydoc alphabet_tuple_base::operator=(component_type const alph)
    SEQAN3_DOXYGEN_ONLY((constexpr qualified & operator=(component_type const alph) noexcept {}))
    //!\copydoc alphabet_tuple_base::operator=(indirect_component_type const alph)
    SEQAN3_DOXYGEN_ONLY((constexpr qualified & operator=(indirect_component_type const alph) noexcept {}))
    //!\}

    /*!\name Write functions
     * \{
     */
    /*!\brief Assign from a character. This modifies the internal sequence letter.
     * \details
     * \stableapi{Since version 3.1.}
     */
    constexpr qualified & assign_char(char_type const c) noexcept
    {
        base_type::assign_rank((seqan3::to_rank(seqan3::assign_char_to(c, sequence_alphabet_type{}))
                                * base_type::cummulative_alph_sizes[0])
                               + (base_type::template to_component_rank<1>() * base_type::cummulative_alph_sizes[1]));

        // The above is noticeably faster than (no subtraction and no division):
        // base_type::template assign_component_rank<0>(
        //     seqan3::to_rank(seqan3::assign_char_to(c, sequence_alphabet_type{})));
        return *this;
    }

    /*!\brief Assign from a Phred score value. This modifies the internal quality letter.
     * \details
     * \stableapi{Since version 3.1.}
     */
    constexpr qualified & assign_phred(phred_type const c) noexcept
    {
        seqan3::assign_phred_to(c, get<1>(*this));
        return *this;
    }
    //!\}

    /*!\name Read functions
     * \{
     */
    /*!\brief Return the Phred score value. This reads the internal quality letter.
     * \details
     * \stableapi{Since version 3.1.}
     */
    constexpr phred_type to_phred() const noexcept
    {
        return rank_to_phred[to_rank()];
    }

    /*!\brief Return a character. This reads the internal sequence letter.
     * \details
     * \stableapi{Since version 3.1.}
     */
    constexpr char_type to_char() const noexcept
    {
        return rank_to_char(to_rank());
    }

    /*!\brief Return a qualified where the quality is preserved, but the sequence letter is complemented.
     * \sa seqan3::complement
     * \sa seqan3::nucleotide
     * \details
     * \stableapi{Since version 3.1.}
     */
    constexpr qualified complement() const noexcept
        requires nucleotide_alphabet<sequence_alphabet_t>
    {
        return qualified{seqan3::complement(get<0>(*this)), get<1>(*this)};
    }
    //!\}

    /*!\brief Validate whether a character is valid in the sequence alphabet.
     * \details
     * \experimentalapi{Experimental since version 3.1.}
     */
    static constexpr bool char_is_valid(char_type const c) noexcept
    {
        return char_is_valid_for<sequence_alphabet_type>(c);
    }

private:
    //!\copydoc seqan3::dna4::rank_to_char
    static constexpr char_type rank_to_char(typename base_type::rank_type const rank)
    {
        return rank_to_char_table[rank];
    }

    //!\copydoc seqan3::dna4::rank_to_char_table
    static constexpr std::array<char_type, alphabet_size> rank_to_char_table{
        []() constexpr
        {
            std::array<char_type, alphabet_size> ret{};

            for (size_t i = 0; i < alphabet_size; ++i)
            {
                size_t const seq_rank =
                    (i / base_type::cummulative_alph_sizes[0]) % seqan3::alphabet_size<quality_alphabet_type>;

                ret[i] = seqan3::to_char(seqan3::assign_rank_to(seq_rank, sequence_alphabet_type{}));
            }

            return ret;
        }()};

    //!\brief Rank to Phred conversion table.
    static constexpr std::array<char_type, alphabet_size> rank_to_phred{
        []() constexpr
        {
            std::array<char_type, alphabet_size> ret{};

            for (size_t i = 0; i < alphabet_size; ++i)
            {
                size_t qual_rank =
                    (i / base_type::cummulative_alph_sizes[1]) % seqan3::alphabet_size<quality_alphabet_type>;

                ret[i] = seqan3::to_phred(seqan3::assign_rank_to(qual_rank, quality_alphabet_type{}));
            }

            return ret;
        }()};
};

/*!\brief Type deduction guide enables usage of qualified without specifying template args.
 * \relates qualified
 */
template <typename sequence_alphabet_type, typename quality_alphabet_type>
qualified(sequence_alphabet_type &&, quality_alphabet_type &&)
    -> qualified<std::decay_t<sequence_alphabet_type>, std::decay_t<quality_alphabet_type>>;

} // namespace seqan3
