// -----------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2020, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2020, Knut Reinert & MPI für molekulare Genetik
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at: https://github.com/seqan/seqan3/blob/master/LICENSE.md
// -----------------------------------------------------------------------------------------------------

/*!\file
 * \author Hannes Hauswedell <hannes.hauswedell AT fu-berlin.de>
 * \brief Provides seqan3::rna4, container aliases and string literals.
 */

#pragma once

#include <vector>

#include <seqan3/alphabet/nucleotide/dna4.hpp>
#include <seqan3/alphabet/nucleotide/nucleotide_base.hpp>

// ------------------------------------------------------------------
// rna4
// ------------------------------------------------------------------

namespace seqan3
{

/*!\brief The four letter RNA alphabet of A,C,G,U.
 * \ingroup nucleotide
 * \implements seqan3::nucleotide_alphabet
 * \implements seqan3::writable_alphabet
 * \if DEV \implements seqan3::detail::writable_constexpr_alphabet \endif
 * \implements seqan3::trivially_copyable
 * \implements seqan3::standard_layout
 * \implements std::regular
 *
 * \details
 * This alphabet has the same internal representation as seqan3::dna4, the only difference is that it prints 'U' on
 * character conversion instead of 'T'. You can assign between values of seqan3::dna4 and seqan3::rna4.
 *
 * Like most alphabets, this alphabet cannot be initialised directly from its character representation.
 * Instead initialise/assign from the character literal or use the
 * function seqan3::rna4::assign_char().
 *
 *\include test/snippet/alphabet/nucleotide/rna4.cpp
 */
class rna4 : public nucleotide_base<rna4, 4>
{
private:
    //!\brief The base class.
    using base_t = nucleotide_base<rna4, 4>;

    //!\brief Befriend seqan3::nucleotide_base.
    friend base_t;
    //!\cond \brief Befriend seqan3::alphabet_base.
    friend base_t::base_t;
    //!\endcond

public:
    /*!\name Constructors, destructor and assignment
     * \{
     */
    constexpr rna4()                          noexcept = default; //!< Defaulted.
    constexpr rna4(rna4 const &)              noexcept = default; //!< Defaulted.
    constexpr rna4(rna4 &&)                   noexcept = default; //!< Defaulted.
    constexpr rna4 & operator=(rna4 const &)  noexcept = default; //!< Defaulted.
    constexpr rna4 & operator=(rna4 &&)       noexcept = default; //!< Defaulted.
    ~rna4()                                   noexcept = default; //!< Defaulted.

    using base_t::base_t;

    //!\brief Allow implicit construction from dna/rna of the same size.
    constexpr rna4(dna4 const & r) noexcept
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
        'C',
        'G',
        'U'
    };

    //!\copydoc seqan3::dna4::char_to_rank
    static constexpr std::array<rank_type, 256> char_to_rank = dna4::char_to_rank;

    //!\copydoc seqan3::dna4::complement_table
    static const std::array<rna4, alphabet_size> complement_table;
};

// ------------------------------------------------------------------
// containers
// ------------------------------------------------------------------

//!\brief Alias for an std::vector of seqan3::rna4.
//!\relates rna4
using rna4_vector = std::vector<rna4>;

// ------------------------------------------------------------------
// literals
// ------------------------------------------------------------------

/*!\name Literals
 * \{
 */

/*!\brief The seqan3::rna4 char literal.
 * \relates seqan3::rna4
 * \returns seqan3::rna4
 */
constexpr rna4 operator""_rna4(char const c) noexcept
{
    return rna4{}.assign_char(c);
}

/*!\brief The seqan3::rna4 string literal.
 * \relates seqan3::rna4
 * \returns seqan3::rna4_vector
 *
 * You can use this string literal to easily assign to rna4_vector:
 *
 * \include test/snippet/alphabet/nucleotide/rna4_literal.cpp
 *
 */
inline rna4_vector operator""_rna4(char const * s, std::size_t n)
{
    rna4_vector r;
    r.resize(n);

    for (size_t i = 0; i < n; ++i)
        r[i].assign_char(s[i]);

    return r;
}
//!\}

// ------------------------------------------------------------------
// rna4 (deferred definition)
// ------------------------------------------------------------------

constexpr std::array<rna4, rna4::alphabet_size> rna4::complement_table
{
    'U'_rna4,    // complement of 'A'_rna4
    'G'_rna4,    // complement of 'C'_rna4
    'C'_rna4,    // complement of 'G'_rna4
    'A'_rna4     // complement of 'U'_rna4
};

} // namespace seqan3
