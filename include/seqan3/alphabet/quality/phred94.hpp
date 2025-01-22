// SPDX-FileCopyrightText: 2006-2025 Knut Reinert & Freie Universität Berlin
// SPDX-FileCopyrightText: 2016-2025 Knut Reinert & MPI für molekulare Genetik
// SPDX-License-Identifier: BSD-3-Clause

/*!\file
 * \author Lydia Buntrock <lydia.buntrock AT fu-berlin.de>
 * \brief Provides seqan3::phred94 quality scores.
 */

#pragma once

#include <seqan3/alphabet/quality/phred_base.hpp>

namespace seqan3
{

/*!\brief Quality type for PacBio Phred scores of HiFi reads.
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
 * The phred94 quality alphabet represents the zero-based Phred score range [0..93] mapped to the ASCII range
 * ['!' .. '~'] (Sanger, Illumina 1.8+ format). It is typically used for HiFi reads produced by PacBio. For Sanger and
 * Illumina Phred scores of raw reads the range is typically (0 to 41), represented as seqan3::phred42. If you expect
 * only slightly larger score types you can use seqan3::phred63 (0 to 62) which still has memory advantages when used
 * with seqan3::qualified.
 *
 * \include test/snippet/alphabet/quality/phred94.cpp
 *
 * \see quality, it contains a comprehensive overview over all Phred implementations and what a Phred score represents.
 *
 * \stableapi{Since version 3.1.}
 */
class phred94 : public phred_base<phred94, 94>
{
private:
    //!\brief The base class.
    using base_t = phred_base<phred94, 94>;

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
    constexpr phred94() noexcept = default;                            //!< Defaulted.
    constexpr phred94(phred94 const &) noexcept = default;             //!< Defaulted.
    constexpr phred94(phred94 &&) noexcept = default;                  //!< Defaulted.
    constexpr phred94 & operator=(phred94 const &) noexcept = default; //!< Defaulted.
    constexpr phred94 & operator=(phred94 &&) noexcept = default;      //!< Defaulted.
    ~phred94() noexcept = default;                                     //!< Defaulted.

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
/*!\brief The seqan3::phred94 char literal.
 * \relatesalso seqan3::phred94
 * \returns seqan3::phred94
 * \details
 *
 * You can use this char literal to assign a seqan3::phred94 character:
 * \include test/snippet/alphabet/quality/phred94_char_literal.cpp
 *
 * \stableapi{Since version 3.1.}
 */
constexpr phred94 operator""_phred94(char const c) noexcept
{
    return phred94{}.assign_char(c);
}

/*!\brief The seqan3::phred94 string literal.
 * \param[in] s A pointer to the character sequence to assign from.
 * \param[in] n The length of the character sequence to assign from.
 * \relates seqan3::phred94
 * \returns seqan3::std::vector<seqan3::phred94>
 *
 * You can use this string literal to easily assign to std::vector<seqan3::phred94>:
 * \include test/snippet/alphabet/quality/phred94_literal.cpp
 *
 * \stableapi{Since version 3.1.}
 */
SEQAN3_WORKAROUND_LITERAL std::vector<phred94> operator""_phred94(char const * s, std::size_t n)
{
    std::vector<phred94> r;
    r.resize(n);

    for (size_t i = 0; i < n; ++i)
        r[i].assign_char(s[i]);

    return r;
}
//!\}

} // namespace literals

} // namespace seqan3
