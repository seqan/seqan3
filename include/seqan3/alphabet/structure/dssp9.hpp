// -----------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2020, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2020, Knut Reinert & MPI für molekulare Genetik
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at: https://github.com/seqan/seqan3/blob/master/LICENSE.md
// -----------------------------------------------------------------------------------------------------

/*!\file
 * \author Joerg Winkler <j.winkler AT fu-berlin.de>
 * \brief Provides the dssp format for protein structure.
 */

#pragma once

#include <vector>

#include <seqan3/alphabet/concept.hpp>
#include <seqan3/alphabet/alphabet_base.hpp>
#include <seqan3/core/char_operations/transform.hpp>

// ------------------------------------------------------------------
// dssp9
// ------------------------------------------------------------------

namespace seqan3
{

/*!\brief The protein structure alphabet of the characters "HGIEBTSCX".
 * \implements seqan3::writable_alphabet
 * \implements seqan3::writable_alphabet
 * \if DEV \implements seqan3::detail::writable_constexpr_alphabet \endif
 * \implements seqan3::trivially_copyable
 * \implements seqan3::standard_layout
 * \implements std::regular
 *
 * \ingroup structure
 *
 * \details
 * The DSSP annotation links structure elements to protein sequences.
 * Originally created with 7 letters as a file format for the DSSP program (http://www.cmbi.ru.nl/dssp.html),
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
 * ###Usage
 * The following code example creates a dssp9 vector, modifies it, and prints the result to stderr.
 * \include test/snippet/alphabet/structure/dssp9_general.cpp
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
    constexpr dssp9()                           noexcept = default; //!< Defaulted.
    constexpr dssp9(dssp9 const &)              noexcept = default; //!< Defaulted.
    constexpr dssp9(dssp9 &&)                   noexcept = default; //!< Defaulted.
    constexpr dssp9 & operator=(dssp9 const &)  noexcept = default; //!< Defaulted.
    constexpr dssp9 & operator=(dssp9 &&)       noexcept = default; //!< Defaulted.
    ~dssp9()                                    noexcept = default; //!< Defaulted.
    //!\}

protected:
    //!\privatesection

    //!\brief Value-to-char conversion table.
    static constexpr char_type rank_to_char[alphabet_size]
    {
        'H', 'B', 'E', 'G', 'I', 'T', 'S', 'C', 'X'
    };

    //!\brief Char-to-value conversion table.
    static constexpr std::array<rank_type, 256> char_to_rank
    {
        [] () constexpr
        {
            std::array<rank_type, 256> ret{};

            // initialize with X (std::array::fill unfortunately not constexpr)
            for (rank_type & rnk : ret)
                rnk = 8u;

            // reverse mapping for characters
            for (rank_type rnk = 0u; rnk < alphabet_size; ++rnk)
            {
                ret[static_cast<rank_type>(rank_to_char[rnk])] = rnk;
            }

            return ret;
        } ()
    };
};

/*!\name Literals
 * \{
 */
/*!\brief The seqan3::dssp9 string literal.
 * \relates seqan3::dssp9
 * \param[in] str A pointer to the character string to assign.
 * \param[in] len The size of the character string to assign.
 * \returns std::vector<seqan3::dssp9>
 *
 * You can use this string literal to easily assign to a vector of seqan3::dssp9 characters:
 * \include test/snippet/alphabet/structure/dssp9_literal.cpp
 */
inline std::vector<dssp9> operator""_dssp9(const char * str, std::size_t len)
{
    std::vector<dssp9> vec;
    vec.resize(len);

    for (size_t idx = 0u; idx < len; ++idx)
        vec[idx].assign_char(str[idx]);

    return vec;
}

/*!\brief The seqan3::dssp9 char literal.
 * \relates seqan3::dssp9
 * \param[in] ch The character to represent as dssp.
 * \returns seqan3::dssp9
 *
 * You can use this string literal to assign a seqan3::dssp9 character:
 * \include test/snippet/alphabet/structure/dssp9_char_literal.cpp
 */
constexpr dssp9 operator""_dssp9(char const ch) noexcept
{
    return dssp9{}.assign_char(ch);
}

//!\}

} // namespace seqan3
