// -----------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2020, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2020, Knut Reinert & MPI für molekulare Genetik
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at: https://github.com/seqan/seqan3/blob/master/LICENSE.md
// -----------------------------------------------------------------------------------------------------

/*!\file
 * \author Svenja Mehringer <svenja.mehringer AT fu-berlin.de>
 * \brief Provides seqan3::dna16sam.
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
 * The seqan3::dna16sam alphabet is the nucleotide alphabet used inside the SAM, BAM and CRAM formats.
 * It has all the letters of the seqan3::dna15 alphabet and the extra alphabet character '=' which denotes a
 * nucleotide character identical to the reference.
 * Without the context of this reference sequence, no assumptions can be made about the actual value of '=' letter.
 *
 * Note that you can assign 'U' as a character to dna16sam and it will silently
 * be converted to 'T'.
 * Lower case letters are accepted when assigning from char (just like seqan3::dna15) and unknown characters are
 * silently converted to 'N'.
 *
 * The complement is the same as for seqan3::dna15, with the addition that the complement of '=' is unknown and
 * therefore set to 'N'.
 *
 * \include test/snippet/alphabet/nucleotide/dna16sam.cpp
 *
 * \stableapi{Since version 3.1.}
 */
class dna16sam : public nucleotide_base<dna16sam, 16>
{
private:
    //!\brief The base class.
    using base_t = nucleotide_base<dna16sam, 16>;

    //!\brief Befriend seqan3::nucleotide_base.
    friend base_t;
    //!\cond \brief Befriend seqan3::alphabet_base.
    friend base_t::base_t;
    //!\endcond

public:
    /*!\name Constructors, destructor and assignment
     * \{
     */
    constexpr dna16sam()                             noexcept = default; //!< Defaulted.
    constexpr dna16sam(dna16sam const &)             noexcept = default; //!< Defaulted.
    constexpr dna16sam(dna16sam &&)                  noexcept = default; //!< Defaulted.
    constexpr dna16sam & operator=(dna16sam const &) noexcept = default; //!< Defaulted.
    constexpr dna16sam & operator=(dna16sam &&)      noexcept = default; //!< Defaulted.
    ~dna16sam()                                      noexcept = default; //!< Defaulted.

    using base_t::base_t;
    //!\}

private:
    //!\copydoc seqan3::dna4::rank_to_char_table
    static constexpr char_type rank_to_char_table[alphabet_size]
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

    //!\copydoc seqan3::dna4::char_to_rank_table
    static constexpr std::array<rank_type, 256> char_to_rank_table
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
                ret[rank_to_char_table[rnk]] = rnk;
                ret[to_lower(rank_to_char_table[rnk])] = rnk;
            }

            // set U equal to T
            ret['U'] = ret['T']; ret['u'] = ret['t'];

            return ret;
        }()
    };

    //!\copydoc seqan3::dna4::complement_table
    static const std::array<dna16sam, alphabet_size> complement_table;

    /*!\copydoc seqan3::dna4::rank_to_char
     *
     * The representation is the same as in the SAM specifications (which is NOT in alphabetical order).
     */
    static constexpr char_type rank_to_char(rank_type const rank)
    {
        return rank_to_char_table[rank];
    }

    //!\copydoc seqan3::dna4::char_to_rank
    static constexpr rank_type char_to_rank(char_type const chr)
    {
        using index_t = std::make_unsigned_t<char_type>;
        return char_to_rank_table[static_cast<index_t>(chr)];
    }
};

// ------------------------------------------------------------------
// containers
// ------------------------------------------------------------------

/*!\brief Alias for an std::vector of seqan3::dna16sam.
 * \relates dna16sam
 * \details
 * \stableapi{Since version 3.1.}
 */
using dna16sam_vector = std::vector<dna16sam>;

// ------------------------------------------------------------------
// literals
// ------------------------------------------------------------------

/*!\name Literals
 * \{
 */

/*!\brief The seqan3::dna16sam char literal.
 * \relates seqan3::dna16sam
 * \returns seqan3::dna16sam
 * \param[in] c The character to assign from.
 * \details
 * \stableapi{Since version 3.1.}
 */
constexpr dna16sam operator""_dna16sam(char const c) noexcept
{
    return dna16sam{}.assign_char(c);
}

/*!\brief The seqan3::dna16sam string literal.
 * \relates seqan3::dna16sam
 * \returns seqan3::dna16sam_vector
 * \param[in] s The string literal to assign from.
 * \param[in] n The length of the string literal s.
 *
 * You can use this string literal to easily assign to seqan3::dna16sam_vector:
 *
 * \include test/snippet/alphabet/nucleotide/dna16sam_literal.cpp
 *
 * \stableapi{Since version 3.1.}
 */
inline dna16sam_vector operator""_dna16sam(char const * s, size_t n)
{
    dna16sam_vector r;
    r.resize(n);

    for (size_t i = 0; i < n; ++i)
        r[i].assign_char(s[i]);

    return r;
}
//!\}

// ------------------------------------------------------------------
// complement deferred definition
// ------------------------------------------------------------------

constexpr std::array<dna16sam, dna16sam::alphabet_size> dna16sam::complement_table
{
    'N'_dna16sam,   // complement of '='_dna16sam
    'T'_dna16sam,   // complement of 'A'_dna16sam
    'G'_dna16sam,   // complement of 'C'_dna16sam
    'K'_dna16sam,   // complement of 'M'_dna16sam
    'C'_dna16sam,   // complement of 'G'_dna16sam
    'Y'_dna16sam,   // complement of 'R'_dna16sam
    'S'_dna16sam,   // complement of 'S'_dna16sam
    'B'_dna16sam,   // complement of 'V'_dna16sam
    'A'_dna16sam,   // complement of 'T'_dna16sam
    'W'_dna16sam,   // complement of 'W'_dna16sam
    'R'_dna16sam,   // complement of 'Y'_dna16sam
    'D'_dna16sam,   // complement of 'H'_dna16sam
    'M'_dna16sam,   // complement of 'K'_dna16sam
    'H'_dna16sam,   // complement of 'D'_dna16sam
    'V'_dna16sam,   // complement of 'B'_dna16sam
    'N'_dna16sam    // complement of 'N'_dna16sam
};

} // namespace seqan3

namespace seqan3
{
//!\deprecated Please use seqan3::dna16sam instead.
using sam_dna16 SEQAN3_DEPRECATED_310 = seqan3::dna16sam;

//!\deprecated Please use seqan3::dna16sam_vector instead.
using sam_dna16_vector SEQAN3_DEPRECATED_310 = dna16sam_vector;

//!\deprecated Please use seqan3::operator""_dna16sam instead.
SEQAN3_DEPRECATED_310 constexpr dna16sam operator""_sam_dna16(char const c) noexcept
{
    return seqan3::operator""_dna16sam(c);
}

//!\deprecated Please use seqan3::operator""_dna16sam instead.
SEQAN3_DEPRECATED_310 inline dna16sam_vector operator""_sam_dna16(char const * s, size_t n)
{
    return seqan3::operator""_dna16sam(s, n);
}
} // namespace seqan3
