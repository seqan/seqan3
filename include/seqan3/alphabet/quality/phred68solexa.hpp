// SPDX-FileCopyrightText: 2006-2025 Knut Reinert & Freie Universität Berlin
// SPDX-FileCopyrightText: 2016-2025 Knut Reinert & MPI für molekulare Genetik
// SPDX-License-Identifier: BSD-3-Clause

/*!\file
 * \author Marie Hoffmann <marie.hoffmann AT fu-berlin.de>
 * \brief Provides seqan3::phred68solexa quality scores.
 */

#pragma once

#include <seqan3/alphabet/quality/phred_base.hpp>

namespace seqan3
{

/*!\brief Quality type for Solexa and deprecated Illumina formats.
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
 * The phred68solexa quality alphabet represents the -5-based Phred score range
 * [-5..62] mapped to the ASCII range [';' .. '~']. It represents the Solexa and
 * the Illumina [1.0;1.8[ standard.
 *
 * \include test/snippet/alphabet/quality/phred68solexa.cpp
 *
 * \stableapi{Since version 3.1.}
 */
class phred68solexa : public phred_base<phred68solexa, 68>
{
private:
    //!\brief The base class.
    using base_t = phred_base<phred68solexa, 68>;

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
    constexpr phred68solexa() noexcept = default;                                  //!< Defaulted.
    constexpr phred68solexa(phred68solexa const &) noexcept = default;             //!< Defaulted.
    constexpr phred68solexa(phred68solexa &&) noexcept = default;                  //!< Defaulted.
    constexpr phred68solexa & operator=(phred68solexa const &) noexcept = default; //!< Defaulted.
    constexpr phred68solexa & operator=(phred68solexa &&) noexcept = default;      //!< Defaulted.
    ~phred68solexa() noexcept = default;                                           //!< Defaulted.

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
    static constexpr phred_type offset_phred{-5};

    /*!\brief The projection offset between char and rank score representation.
     * \details
     * \stableapi{Since version 3.1.}
     */
    static constexpr char_type offset_char{';'};
    //!\}
};

inline namespace literals
{

/*!\name Quality literals
 * \{
 */
/*!\brief The seqan3::phred68solexa char literal.
 * \relatesalso seqan3::phred68solexa
 * \returns seqan3::phred68solexa
 * \details
 *
 * You can use this char literal to assign a seqan3::phred68solexa character:
 * \include test/snippet/alphabet/quality/phred68solexa_char_literal.cpp
 *
 * \stableapi{Since version 3.1.}
 */
constexpr phred68solexa operator""_phred68solexa(char const c) noexcept
{
    return phred68solexa{}.assign_char(c);
}

/*!\brief The seqan3::phred68solexa string literal.
 * \param[in] s A pointer to the character sequence to assign from.
 * \param[in] n The length of the character sequence to assign from.
 * \relates seqan3::phred68solexa
 * \returns seqan3::std::vector<seqan3::phred68solexa>
 *
 * You can use this string literal to easily assign to std::vector<seqan3::phred68solexa>:
 * \include test/snippet/alphabet/quality/phred68solexa_literal.cpp
 *
 * \stableapi{Since version 3.1.}
 */
SEQAN3_WORKAROUND_LITERAL std::vector<phred68solexa> operator""_phred68solexa(char const * s, std::size_t n)
{
    std::vector<phred68solexa> r;
    r.resize(n);

    for (size_t i = 0; i < n; ++i)
        r[i].assign_char(s[i]);

    return r;
}
//!\}

} // namespace literals

} // namespace seqan3
