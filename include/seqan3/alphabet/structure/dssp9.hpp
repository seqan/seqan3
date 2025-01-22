// SPDX-FileCopyrightText: 2006-2025 Knut Reinert & Freie Universität Berlin
// SPDX-FileCopyrightText: 2016-2025 Knut Reinert & MPI für molekulare Genetik
// SPDX-License-Identifier: BSD-3-Clause

/*!\file
 * \author Joerg Winkler <j.winkler AT fu-berlin.de>
 * \brief Provides the dssp format for protein structure.
 */

#pragma once

#include <vector>

#include <seqan3/alphabet/alphabet_base.hpp>
#include <seqan3/alphabet/concept.hpp>
#include <seqan3/utility/char_operations/transform.hpp>

// ------------------------------------------------------------------
// dssp9
// ------------------------------------------------------------------

namespace seqan3
{

/*!\brief The protein structure alphabet of the characters "HGIEBTSCX".
 * \implements seqan3::writable_alphabet
 * \implements seqan3::detail::writable_constexpr_alphabet
 * \implements seqan3::trivially_copyable
 * \implements seqan3::standard_layout
 * \implements std::regular
 *
 * \ingroup alphabet_structure
 *
 * \details
 * The DSSP annotation links structure elements to protein sequences.
 * Originally created with 7 letters as a file format for the DSSP program (https://swift.cmbi.umcn.nl/gv/dssp/),
 * it is also used in the stockholm file format for structure alignments, extended by the characters C and X
 * (https://en.wikipedia.org/wiki/Stockholm_format).
 *
 * The letter abbreviations are as follows:
 *
 * H = alpha helix
 * B = beta bridge
 * E = strand
 * G = helix-3
 * I = helix-5
 * T = turn
 * S = bend
 * C = coil/loop
 * X = unknown
 *
 * ### Example
 *
 * \include test/snippet/alphabet/structure/dssp9.cpp
 *
 * \experimentalapi{Experimental since version 3.1.}
 */
class dssp9 : public alphabet_base<dssp9, 9>
{
private:
    //!\brief The base class.
    using base_t = alphabet_base<dssp9, 9>;

    //!\brief Befriend seqan3::alphabet_base.
    friend base_t;

public:
    /*!\name Constructors, destructor and assignment
     * \{
     */
    constexpr dssp9() noexcept = default;                          //!< Defaulted.
    constexpr dssp9(dssp9 const &) noexcept = default;             //!< Defaulted.
    constexpr dssp9(dssp9 &&) noexcept = default;                  //!< Defaulted.
    constexpr dssp9 & operator=(dssp9 const &) noexcept = default; //!< Defaulted.
    constexpr dssp9 & operator=(dssp9 &&) noexcept = default;      //!< Defaulted.
    ~dssp9() noexcept = default;                                   //!< Defaulted.

    //!\}

private:
    //!\copydoc seqan3::dna4::rank_to_char_table
    static constexpr char_type rank_to_char_table[alphabet_size]{'H', 'B', 'E', 'G', 'I', 'T', 'S', 'C', 'X'};

    //!\copydoc seqan3::dna4::rank_to_char
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
    static constexpr std::array<rank_type, 256> char_to_rank_table{
        []() constexpr
        {
            std::array<rank_type, 256> ret{};

            ret.fill(8u);

            // reverse mapping for characters
            for (rank_type rnk = 0u; rnk < alphabet_size; ++rnk)
                ret[static_cast<rank_type>(rank_to_char_table[rnk])] = rnk;

            return ret;
        }()};
};

inline namespace literals
{

/*!\name Structure literals
 * \{
 */
/*!\brief The seqan3::dssp9 char literal.
 * \relatesalso seqan3::dssp9
 * \param[in] ch The character to represent as dssp.
 * \returns seqan3::dssp9
 *
 * You can use this char literal to assign a seqan3::dssp9 character:
 * \include test/snippet/alphabet/structure/dssp9_char_literal.cpp
 *
 * \experimentalapi{Experimental since version 3.1.}
 */
constexpr dssp9 operator""_dssp9(char const ch) noexcept
{
    return dssp9{}.assign_char(ch);
}

/*!\brief The seqan3::dssp9 string literal.
 * \relatesalso seqan3::dssp9
 * \param[in] str A pointer to the character string to assign.
 * \param[in] len The size of the character string to assign.
 * \returns std::vector<seqan3::dssp9>
 *
 * You can use this string literal to easily assign to std::vector<seqan3::dssp9>:
 * \include test/snippet/alphabet/structure/dssp9_literal.cpp
 *
 * \experimentalapi{Experimental since version 3.1.}
 */
SEQAN3_WORKAROUND_LITERAL std::vector<dssp9> operator""_dssp9(char const * str, std::size_t len)
{
    std::vector<dssp9> vec;
    vec.resize(len);

    for (size_t idx = 0u; idx < len; ++idx)
        vec[idx].assign_char(str[idx]);

    return vec;
}
//!\}

} // namespace literals

} // namespace seqan3
