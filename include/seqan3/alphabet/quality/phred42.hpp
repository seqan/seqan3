// SPDX-FileCopyrightText: 2006-2025 Knut Reinert & Freie Universität Berlin
// SPDX-FileCopyrightText: 2016-2025 Knut Reinert & MPI für molekulare Genetik
// SPDX-License-Identifier: BSD-3-Clause

/*!\file
 * \author Marie Hoffmann <marie.hoffmann AT fu-berlin.de>
 * \brief Provides seqan3::phred42 quality scores.
 */

#pragma once

#include <seqan3/alphabet/quality/phred_base.hpp>

namespace seqan3
{

/*!\brief Quality type for traditional Sanger and modern Illumina Phred scores.
 * \implements seqan3::writable_quality_alphabet
 * \implements seqan3::detail::writable_constexpr_alphabet
 * \implements seqan3::trivially_copyable
 * \implements seqan3::standard_layout
 * \implements std::regular
 *
 * \ingroup alphabet_quality
 *
 * \details
 *
 * The phred42 \ref alphabet_quality alphabet represents the zero-based Phred score range [0..41] mapped to the
 * consecutive ASCII range ['!' .. 'J']. It therefore can represent the Illumina 1.8+ standard and the original Sanger
 * score. If you intend to use Phred scores exceeding 41, use the larger score types, namely seqan3::phred63 or
 * seqan3::phred94, otherwise on construction exceeding scores are mapped to 41.
 * Via seqan3::qualified, you can combine a nucleotide alphabet with the Phred score to save space.
 * All dna and rna combinations with seqan3::phred42 still fit into a single byte,
 * e.g. `seqan3::qualified<seqan3::dna4, seqan3::phred42>` (4 * 42 = 168 values can be stored in a single byte which can
 * contain up to 256 values).
 *
 * \include test/snippet/alphabet/quality/phred42.cpp
 *
 * \see quality, it contains a comprehensive overview over all Phred implementations and what a Phred score represents.
 *
 * \stableapi{Since version 3.1.}
 */
class phred42 : public phred_base<phred42, 42>
{
private:
    //!\brief The base class.
    using base_t = phred_base<phred42, 42>;

    //!\brief Befriend seqan3::phred_base.
    friend base_t;
    //!\cond
    //!\brief Befriend seqan3::alphabet_base.
    friend base_t::base_t;
    //!\endcond

public:
    /*!\name Constructors, destructor and assignment
     * \{
     */
    constexpr phred42() noexcept = default;                            //!< Defaulted.
    constexpr phred42(phred42 const &) noexcept = default;             //!< Defaulted.
    constexpr phred42(phred42 &&) noexcept = default;                  //!< Defaulted.
    constexpr phred42 & operator=(phred42 const &) noexcept = default; //!< Defaulted.
    constexpr phred42 & operator=(phred42 &&) noexcept = default;      //!< Defaulted.
    ~phred42() noexcept = default;                                     //!< Defaulted.

    // Inherit converting constructor
    using base_t::base_t;
    //!\}

    /*!\name Member variables.
     * \{
     */
    /*!\brief The projection offset between Phred and rank score representation.
     * \details
     * \stableapi{Since version 3.1.}
     */
    static constexpr phred_type offset_phred{0};

    /*!\brief The projection offset between char and rank score representation.
     * \details
     * \stableapi{Since version 3.1.}
     */
    static constexpr char_type offset_char{'!'};
    //!\}
};

inline namespace literals
{

/*!\name Quality literals
 * \{
 */
/*!\brief The seqan3::phred42 char literal.
 * \relatesalso seqan3::phred42
 * \returns seqan3::phred42
 * \details
 *
 * You can use this char literal to assign a seqan3::phred42 character:
 * \include test/snippet/alphabet/quality/phred42_char_literal.cpp
 *
 * \stableapi{Since version 3.1.}
 */
constexpr phred42 operator""_phred42(char const c) noexcept
{
    return phred42{}.assign_char(c);
}

/*!\brief The seqan3::phred42 string literal.
 * \param[in] s A pointer to the character sequence to assign from.
 * \param[in] n The length of the character sequence to assign from.
 * \relates seqan3::phred42
 * \returns seqan3::std::vector<seqan3::phred42>
 *
 * You can use this string literal to easily assign to std::vector<seqan3::phred42>:
 * \include test/snippet/alphabet/quality/phred42_literal.cpp
 *
 * \stableapi{Since version 3.1.}
 */
SEQAN3_WORKAROUND_LITERAL std::vector<phred42> operator""_phred42(char const * s, std::size_t n)
{
    std::vector<phred42> r;
    r.resize(n);

    for (size_t i = 0; i < n; ++i)
        r[i].assign_char(s[i]);

    return r;
}
//!\}

} // namespace literals

} // namespace seqan3
