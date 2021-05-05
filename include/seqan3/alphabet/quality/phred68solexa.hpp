// -----------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2020, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2020, Knut Reinert & MPI für molekulare Genetik
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at: https://github.com/seqan/seqan3/blob/master/LICENSE.md
// -----------------------------------------------------------------------------------------------------

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
 * \if DEV \implements seqan3::detail::writable_constexpr_alphabet \endif
 * \implements seqan3::trivially_copyable
 * \implements seqan3::standard_layout
 * \implements std::regular
 *
 * \ingroup quality
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
    //!\cond \brief Befriend seqan3::alphabet_base.
    friend base_t::base_t;
    //!\endcond

public:
    /*!\name Constructors, destructor and assignment
     * \{
     */
    constexpr phred68solexa()                                   noexcept = default; //!< Defaulted.
    constexpr phred68solexa(phred68solexa const &)              noexcept = default; //!< Defaulted.
    constexpr phred68solexa(phred68solexa &&)                   noexcept = default; //!< Defaulted.
    constexpr phred68solexa & operator=(phred68solexa const &)  noexcept = default; //!< Defaulted.
    constexpr phred68solexa & operator=(phred68solexa &&)       noexcept = default; //!< Defaulted.
    ~phred68solexa()                                            noexcept = default; //!< Defaulted.

    /*!\brief Allow construction from the Phred score value.
     * \details
     * \deprecated This will be removed in 3.1.0. Please use seqan3::phred68solexa::assign_phred() or '!'_phred68solexa.
     */
    SEQAN3_DEPRECATED_310 constexpr phred68solexa(phred_type const p)
    {
        assign_phred(p);
    }

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

//!\deprecated Please use seqan3::phred68solexa instead.
using phred68legacy SEQAN3_DEPRECATED_310 = seqan3::phred68solexa;

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
inline std::vector<phred68solexa> operator""_phred68solexa(char const * s, std::size_t n)
{
    std::vector<phred68solexa> r;
    r.resize(n);

    for (size_t i = 0; i < n; ++i)
        r[i].assign_char(s[i]);

    return r;
}

//!\deprecated Please use seqan3::operator""_phred68solexa instead.
SEQAN3_DEPRECATED_310 constexpr phred68solexa operator""_phred68legacy(char const c) noexcept
{
    return seqan3::operator""_phred68solexa(c);
}

//!\deprecated Please use seqan3::operator""_phred68solexa instead.
SEQAN3_DEPRECATED_310 inline std::vector<phred68solexa> operator""_phred68legacy(char const * s, size_t n)
{
    return seqan3::operator""_phred68solexa(s, n);
}
//!\}

} // inline namespace literals

} // namespace seqan3
