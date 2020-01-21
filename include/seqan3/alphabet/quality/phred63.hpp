// -----------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2020, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2020, Knut Reinert & MPI für molekulare Genetik
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at: https://github.com/seqan/seqan3/blob/master/LICENSE.md
// -----------------------------------------------------------------------------------------------------

/*!\file
 * \author Marie Hoffmann <marie.hoffmann AT fu-berlin.de>
 * \brief Provides seqan3::phred63 quality scores.
 */

#pragma once

#include <seqan3/alphabet/quality/quality_base.hpp>

// ------------------------------------------------------------------
// phred63
// ------------------------------------------------------------------

namespace seqan3
{

/*!\brief Quality type for traditional Sanger and modern Illumina Phred scores (full range).
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
 * The phred63 quality alphabet represents the zero-based phred score range
 * [0..62] mapped to the ASCII range ['!' .. '~']. It represents the Sanger and
 * Illumina 1.8+ standard beyond the typical range of raw reads (0 to 41).
 *
 * \include test/snippet/alphabet/quality/phred63.cpp
 */
class phred63 : public quality_base<phred63, 63>
{
private:
    //!\brief The base class.
    using base_t = quality_base<phred63, 63>;

    //!\brief Befriend seqan3::quality_base.
    friend base_t;
    //!\cond \brief Befriend seqan3::alphabet_base.
    friend base_t::base_t;
    //!\endcond

public:
    /*!\name Constructors, destructor and assignment
     * \{
     */
    constexpr phred63()                             noexcept = default; //!< Defaulted.
    constexpr phred63(phred63 const &)              noexcept = default; //!< Defaulted.
    constexpr phred63(phred63 &&)                   noexcept = default; //!< Defaulted.
    constexpr phred63 & operator=(phred63 const &)  noexcept = default; //!< Defaulted.
    constexpr phred63 & operator=(phred63 &&)       noexcept = default; //!< Defaulted.
    ~phred63()                                      noexcept = default; //!< Defaulted.

    //!\brief Construct from phred value.
    constexpr phred63(phred_type const p) : base_t{p} {}

    // Inherit converting constructor
    using base_t::base_t;
    //!\}

    /*!\name Member variables.
     * \{
     */
    //!\brief The projection offset between phred and rank score representation.
    static constexpr phred_type offset_phred{0};

    //!\brief The projection offset between char and rank score representation.
    static constexpr char_type offset_char{'!'};
    //!\}
};

/*!\name Literals
 * \{
 */
/*!\brief The seqan3::phred63 char literal.
 * \relates seqan3::phred63
 * \returns seqan3::phred63
 */
constexpr phred63 operator""_phred63(char const c) noexcept
{
    return phred63{}.assign_char(c);
}

/*!\brief The seqan3::phred63 string literal.
 * \param[in] s A pointer to the character sequence to assign from.
 * \param[in] n The length of the character sequence to assign from.
 * \relates seqan3::phred63
 * \returns seqan3::std::vector<seqan3::phred63>
 *
 * You can use this string literal to easily assign to std::vector<seqan3::phred63>:
 *
 * \include test/snippet/alphabet/quality/phred63_literal.cpp
 */
inline std::vector<phred63> operator""_phred63(char const * s, std::size_t n)
{
    std::vector<phred63> r;
    r.resize(n);

    for (size_t i = 0; i < n; ++i)
        r[i].assign_char(s[i]);

    return r;
}
//!\}

} // namespace seqan3
