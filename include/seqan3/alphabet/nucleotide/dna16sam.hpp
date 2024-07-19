// SPDX-FileCopyrightText: 2006-2024 Knut Reinert & Freie Universität Berlin
// SPDX-FileCopyrightText: 2016-2024 Knut Reinert & MPI für molekulare Genetik
// SPDX-License-Identifier: BSD-3-Clause

/*!\file
 * \author Svenja Mehringer <svenja.mehringer AT fu-berlin.de>
 * \brief Provides seqan3::dna16sam.
 */

#pragma once

#include <seqan3/alphabet/nucleotide/nucleotide_base.hpp>

namespace seqan3
{

/*!\brief A 16 letter DNA alphabet, containing all IUPAC symbols minus the gap and plus an equality sign ('=').
 * \ingroup alphabet_nucleotide
 * \implements seqan3::writable_alphabet
 * \implements seqan3::nucleotide_alphabet
 * \implements seqan3::trivially_copyable
 * \implements seqan3::standard_layout
 * \implements seqan3::detail::writable_constexpr_alphabet
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
    //!\cond
    //!\brief Befriend seqan3::alphabet_base.
    friend base_t::base_t;
    //!\endcond

public:
    /*!\name Constructors, destructor and assignment
     * \{
     */
    constexpr dna16sam() noexcept = default;                             //!< Defaulted.
    constexpr dna16sam(dna16sam const &) noexcept = default;             //!< Defaulted.
    constexpr dna16sam(dna16sam &&) noexcept = default;                  //!< Defaulted.
    constexpr dna16sam & operator=(dna16sam const &) noexcept = default; //!< Defaulted.
    constexpr dna16sam & operator=(dna16sam &&) noexcept = default;      //!< Defaulted.
    ~dna16sam() noexcept = default;                                      //!< Defaulted.

    using base_t::base_t;
    //!\}

private:
    //!\copydoc seqan3::dna4::rank_to_char_table
    static constexpr char_type rank_to_char_table
        [alphabet_size]{'=', 'A', 'C', 'M', 'G', 'R', 'S', 'V', 'T', 'W', 'Y', 'H', 'K', 'D', 'B', 'N'};

    //!\copydoc seqan3::dna4::rank_complement_table
    static constexpr rank_type rank_complement_table[alphabet_size]{
        15, // N is complement of '='_dna16sam  0
        8,  // T is complement of 'A'_dna16sam  1
        4,  // G is complement of 'C'_dna16sam  2
        12, // K is complement of 'M'_dna16sam  3
        2,  // C is complement of 'G'_dna16sam  4
        10, // Y is complement of 'R'_dna16sam  5
        6,  // S is complement of 'S'_dna16sam  6
        14, // B is complement of 'V'_dna16sam  7
        1,  // A is complement of 'T'_dna16sam  8
        9,  // W is complement of 'W'_dna16sam  9
        5,  // R is complement of 'Y'_dna16sam 10
        13, // D is complement of 'H'_dna16sam 11
        3,  // M is complement of 'K'_dna16sam 12
        11, // H is complement of 'D'_dna16sam 13
        7,  // V is complement of 'B'_dna16sam 14
        15  // N is complement of 'N'_dna16sam 15
    };

    //!\copydoc seqan3::dna4::rank_complement
    static constexpr rank_type rank_complement(rank_type const rank)
    {
        return rank_complement_table[rank];
    }

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

    //!\copydoc seqan3::dna4::char_to_rank_table
    static constexpr std::array<rank_type, 256> char_to_rank_table{[]() constexpr
                                                                   {
                                                                       std::array<rank_type, 256> ret{};

                                                                       ret.fill(15u); // initialize with UNKNOWN ('N')

                                                                       // reverse mapping for characters and their lowercase
                                                                       for (size_t rnk = 0u; rnk < alphabet_size; ++rnk)
                                                                       {
                                                                           ret[rank_to_char_table[rnk]] = rnk;
                                                                           ret[to_lower(rank_to_char_table[rnk])] = rnk;
                                                                       }

                                                                       // set U equal to T
                                                                       ret['U'] = ret['T'];
                                                                       ret['u'] = ret['t'];

                                                                       return ret;
                                                                   }()};
};

// ------------------------------------------------------------------
// containers
// ------------------------------------------------------------------

/*!\brief Alias for a std::vector of seqan3::dna16sam.
 * \relates dna16sam
 * \details
 * \stableapi{Since version 3.1.}
 */
using dna16sam_vector = std::vector<dna16sam>;

// ------------------------------------------------------------------
// literals
// ------------------------------------------------------------------
inline namespace literals
{

/*!\name Nucleotide literals
 * \{
 */
/*!\brief The seqan3::dna16sam char literal.
 * \relatesalso seqan3::dna16sam
 * \returns seqan3::dna16sam
 * \param[in] c The character to assign from.
 * \details
 *
 * You can use this char literal to assign a seqan3::dna16sam character:
 * \include test/snippet/alphabet/nucleotide/dna16sam_char_literal.cpp
 *
 * \stableapi{Since version 3.1.}
 */
constexpr dna16sam operator""_dna16sam(char const c) noexcept
{
    return dna16sam{}.assign_char(c);
}

/*!\brief The seqan3::dna16sam string literal.
 * \relatesalso seqan3::dna16sam
 * \returns seqan3::dna16sam_vector
 * \param[in] s The string literal to assign from.
 * \param[in] n The length of the string literal s.
 *
 * You can use this string literal to easily assign to seqan3::dna16sam_vector:
 * \include test/snippet/alphabet/nucleotide/dna16sam_literal.cpp
 *
 * \stableapi{Since version 3.1.}
 */
SEQAN3_WORKAROUND_LITERAL dna16sam_vector operator""_dna16sam(char const * s, size_t n)
{
    dna16sam_vector r;
    r.resize(n);

    for (size_t i = 0; i < n; ++i)
        r[i].assign_char(s[i]);

    return r;
}
//!\}

} // namespace literals

} // namespace seqan3
