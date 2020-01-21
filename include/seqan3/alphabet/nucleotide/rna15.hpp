// -----------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2020, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2020, Knut Reinert & MPI für molekulare Genetik
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at: https://github.com/seqan/seqan3/blob/master/LICENSE.md
// -----------------------------------------------------------------------------------------------------

/*!\file
 * \author Hannes Hauswedell <hannes.hauswedell AT fu-berlin.de>
 * \brief Provides seqan3::rna15, container aliases and string literals.
 */

#pragma once

#include <vector>

#include <seqan3/alphabet/nucleotide/dna15.hpp>
#include <seqan3/alphabet/nucleotide/nucleotide_base.hpp>

// ------------------------------------------------------------------
// rna15
// ------------------------------------------------------------------

namespace seqan3
{

/*!\brief The 15 letter RNA alphabet, containing all IUPAC smybols minus the gap.
 * \implements seqan3::nucleotide_alphabet
 * \implements seqan3::writable_alphabet
 * \if DEV \implements seqan3::detail::writable_constexpr_alphabet \endif
 * \implements seqan3::trivially_copyable
 * \implements seqan3::standard_layout
 * \implements std::regular
 *
 * \ingroup nucleotide
 *
 * \details
 *
 * This alphabet has the same internal representation as seqan3::dna15, the only difference is that it prints 'U' on
 * character conversion instead of 'T'. You can assign between values of seqan3::dna15 and seqan3::rna15.
 *
 * Like most alphabets, this alphabet cannot be initialised directly from its character representation.
 * Instead initialise/assign from the character literal or use the
 * function seqan3::rna15::assign_char().
 *
 *\include test/snippet/alphabet/nucleotide/rna15.cpp
 */
class rna15 : public nucleotide_base<rna15, 15>
{
private:
    //!\brief The base class.
    using base_t = nucleotide_base<rna15, 15>;

    //!\brief Befriend seqan3::nucleotide_base.
    friend base_t;
    //!\cond \brief Befriend seqan3::alphabet_base.
    friend base_t::base_t;
    //!\endcond

public:
    /*!\name Constructors, destructor and assignment
     * \{
     */
    constexpr rna15()                           noexcept = default; //!< Defaulted.
    constexpr rna15(rna15 const &)              noexcept = default; //!< Defaulted.
    constexpr rna15(rna15 &&)                   noexcept = default; //!< Defaulted.
    constexpr rna15 & operator=(rna15 const &)  noexcept = default; //!< Defaulted.
    constexpr rna15 & operator=(rna15 &&)       noexcept = default; //!< Defaulted.
    ~rna15()                                    noexcept = default; //!< Defaulted.

    using base_t::base_t;

    //!\brief Allow implicit construction from dna/rna of the same size.
    constexpr rna15(dna15 const & r) noexcept
#if SEQAN3_WORKAROUND_GCC_90897
        requires true
#endif
    {
        assign_rank(r.to_rank());
    }
    //!\}

protected:
    //!\privatesection

    //!\copydoc seqan3::dna4::rank_to_char
    static constexpr char_type rank_to_char[alphabet_size]
    {
        'A',
        'B',
        'C',
        'D',
        'G',
        'H',
        'K',
        'M',
        'N',
        'R',
        'S',
        'U',
        'V',
        'W',
        'Y'
    };

    //!\copydoc seqan3::dna4::char_to_rank
    static constexpr std::array<rank_type, 256> char_to_rank = dna15::char_to_rank;

    //!\copydoc seqan3::dna4::complement_table
    static const std::array<rna15, alphabet_size> complement_table;
};

// ------------------------------------------------------------------
// containers
// ------------------------------------------------------------------

//!\brief Alias for an std::vector of seqan3::rna15.
//!\relates rna15
using rna15_vector = std::vector<rna15>;

// ------------------------------------------------------------------
// literals
// ------------------------------------------------------------------

/*!\name Literals
 * \{
 */

/*!\brief The seqan3::rna15 char literal.
 * \relates seqan3::rna15
 * \returns seqan3::rna15
 */
constexpr rna15 operator""_rna15(char const c) noexcept
{
    return rna15{}.assign_char(c);
}

/*!\brief The seqan3::rna15 string literal.
 * \relates seqan3::rna15
 * \returns seqan3::rna15_vector
 *
 * You can use this string literal to easily assign to rna15_vector:
 *
 * \include test/snippet/alphabet/nucleotide/rna15_literal.cpp
 *
 */
inline rna15_vector operator""_rna15(char const * s, std::size_t n)
{
    rna15_vector r;
    r.resize(n);

    for (size_t i = 0; i < n; ++i)
        r[i].assign_char(s[i]);

    return r;
}
//!\}

// ------------------------------------------------------------------
// rna15 (deferred definition)
// ------------------------------------------------------------------

constexpr std::array<rna15, rna15::alphabet_size> rna15::complement_table
{
    'U'_rna15,    // complement of 'A'_rna15
    'V'_rna15,    // complement of 'B'_rna15
    'G'_rna15,    // complement of 'C'_rna15
    'H'_rna15,    // complement of 'D'_rna15
    'C'_rna15,    // complement of 'G'_rna15
    'D'_rna15,    // complement of 'H'_rna15
    'M'_rna15,    // complement of 'K'_rna15
    'K'_rna15,    // complement of 'M'_rna15
    'N'_rna15,    // complement of 'N'_rna15
    'Y'_rna15,    // complement of 'R'_rna15
    'S'_rna15,    // complement of 'S'_rna15
    'A'_rna15,    // complement of 'U'_rna15
    'B'_rna15,    // complement of 'V'_rna15
    'W'_rna15,    // complement of 'W'_rna15
    'R'_rna15     // complement of 'Y'_rna15
};

} // namespace seqan3
