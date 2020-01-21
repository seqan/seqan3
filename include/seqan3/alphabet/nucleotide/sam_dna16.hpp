// -----------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2020, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2020, Knut Reinert & MPI für molekulare Genetik
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at: https://github.com/seqan/seqan3/blob/master/LICENSE.md
// -----------------------------------------------------------------------------------------------------

/*!\file
 * \author Svenja Mehringer <svenja.mehringer AT fu-berlin.de>
 * \brief Provides seqan3::sam_dna16.
 */

#pragma once

#include <seqan3/alphabet/nucleotide/nucleotide_base.hpp>

namespace seqan3
{

/*!\brief A 16 letter DNA alphabet, containing all IUPAC symbols minus the gap and plus an equality sign ('=').
 * \ingroup nucleotide
 * \implements seqan3::writable_alphabet
 * \implements seqan3::nucleotide_alphabet
 * \implements seqan3::trivially_copyable
 * \implements seqan3::standard_layout
 * \if DEV \implements seqan3::detail::writable_constexpr_alphabet \endif
 *
 * \details
 *
 * The seqan3::sam_dna16 alphabet is the nucleotide alphabet used inside the SAM, BAM and CRAM formats.
 * It has all the letters of the seqan3::dna15 alphabet and the extra alphabet character '=' which denotes a
 * nucleotide character identical to the reference.
 * Without the context of this reference sequence, no assumptions can be made about the actual value of '=' letter.
 *
 * Note that you can assign 'U' as a character to sam_dna16 and it will silently
 * be converted to 'T'.
 * Lower case letters are accepted when assigning from char (just like seqan3::dna15) and unknown characters are
 * silently converted to 'N'.
 *
 * The complement is the same as for seqan3::dna15, with the addition that the complement of '=' is unknown and
 * therefore set to 'N'.
 *
 * \include test/snippet/alphabet/nucleotide/sam_dna16.cpp
 */
class sam_dna16 : public nucleotide_base<sam_dna16, 16>
{
private:
    //!\brief The base class.
    using base_t = nucleotide_base<sam_dna16, 16>;

    //!\brief Befriend seqan3::nucleotide_base.
    friend base_t;
    //!\cond \brief Befriend seqan3::alphabet_base.
    friend base_t::base_t;
    //!\endcond

public:
    /*!\name Constructors, destructor and assignment
     * \{
     */
    constexpr sam_dna16()                              noexcept = default; //!< Defaulted.
    constexpr sam_dna16(sam_dna16 const &)             noexcept = default; //!< Defaulted.
    constexpr sam_dna16(sam_dna16 &&)                  noexcept = default; //!< Defaulted.
    constexpr sam_dna16 & operator=(sam_dna16 const &) noexcept = default; //!< Defaulted.
    constexpr sam_dna16 & operator=(sam_dna16 &&)      noexcept = default; //!< Defaulted.
    ~sam_dna16()                                       noexcept = default; //!< Defaulted.

    using base_t::base_t;
    //!\}

protected:
    //!\privatesection

    //!\brief The representation is the same as in the SAM specifications (which is NOT in alphabetical order).
    static constexpr char_type rank_to_char[alphabet_size]
    {
        '=',
        'A',
        'C',
        'M',
        'G',
        'R',
        'S',
        'V',
        'T',
        'W',
        'Y',
        'H',
        'K',
        'D',
        'B',
        'N'
    };

    //!\copydoc seqan3::dna4::char_to_rank
    static constexpr std::array<rank_type, 256> char_to_rank
    {
        [] () constexpr
        {
            std::array<rank_type, 256> ret{};

            // initialize with UNKNOWN (std::array::fill unfortunately not constexpr)
            for (auto & c : ret)
                c = 15; // rank of 'N'

            // reverse mapping for characters and their lowercase
            for (size_t rnk = 0u; rnk < alphabet_size; ++rnk)
            {
                ret[         rank_to_char[rnk] ] = rnk;
                ret[to_lower(rank_to_char[rnk])] = rnk;
            }

            // set U equal to T
            ret['U'] = ret['T']; ret['u'] = ret['t'];

            return ret;
        }()
    };

    //!\copydoc seqan3::dna4::complement_table
    static const std::array<sam_dna16, alphabet_size> complement_table;
};

// ------------------------------------------------------------------
// containers
// ------------------------------------------------------------------

//!\brief Alias for an std::vector of seqan3::sam_dna16.
//!\relates sam_dna16
using sam_dna16_vector = std::vector<sam_dna16>;

// ------------------------------------------------------------------
// literals
// ------------------------------------------------------------------

/*!\name Literals
 * \{
 */

/*!\brief The seqan3::sam_dna16 char literal.
 * \relates seqan3::sam_dna16
 * \returns seqan3::sam_dna16
 * \param[in] c The character to assign from.
 */
constexpr sam_dna16 operator""_sam_dna16(char const c) noexcept
{
    return sam_dna16{}.assign_char(c);
}

/*!\brief The seqan3::sam_dna16 string literal.
 * \relates seqan3::sam_dna16
 * \returns seqan3::sam_dna16_vector
 * \param[in] s The string literal to assign from.
 * \param[in] n The length of the string literal s.
 *
 * You can use this string literal to easily assign to seqan3::sam_dna16_vector:
 *
 * \include test/snippet/alphabet/nucleotide/sam_dna16_literal.cpp
 */
inline sam_dna16_vector operator""_sam_dna16(char const * s, size_t n)
{
    sam_dna16_vector r;
    r.resize(n);

    for (size_t i = 0; i < n; ++i)
        r[i].assign_char(s[i]);

    return r;
}
//!\}

// ------------------------------------------------------------------
// complement deferred definition
// ------------------------------------------------------------------

constexpr std::array<sam_dna16, sam_dna16::alphabet_size> sam_dna16::complement_table
{
    'N'_sam_dna16,    // complement of '='_sam_dna16
    'T'_sam_dna16,    // complement of 'A'_sam_dna16
    'G'_sam_dna16,    // complement of 'C'_sam_dna16
    'K'_sam_dna16,    // complement of 'M'_sam_dna16
    'C'_sam_dna16,    // complement of 'G'_sam_dna16
    'Y'_sam_dna16,    // complement of 'R'_sam_dna16
    'S'_sam_dna16,    // complement of 'S'_sam_dna16
    'B'_sam_dna16,    // complement of 'V'_sam_dna16
    'A'_sam_dna16,    // complement of 'T'_sam_dna16
    'W'_sam_dna16,    // complement of 'W'_sam_dna16
    'R'_sam_dna16,    // complement of 'Y'_sam_dna16
    'D'_sam_dna16,    // complement of 'H'_sam_dna16
    'M'_sam_dna16,    // complement of 'K'_sam_dna16
    'H'_sam_dna16,    // complement of 'D'_sam_dna16
    'V'_sam_dna16,    // complement of 'B'_sam_dna16
    'N'_sam_dna16     // complement of 'N'_sam_dna16
};

} // namespace seqan3
