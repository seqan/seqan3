// -----------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2020, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2020, Knut Reinert & MPI für molekulare Genetik
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at: https://github.com/seqan/seqan3/blob/master/LICENSE.md
// -----------------------------------------------------------------------------------------------------

/*!\file
 * \author Hannes Hauswedell <hannes.hauswedell AT fu-berlin.de>
 * \brief Provides seqan3::dna5, container aliases and string literals.
 */

#pragma once

#include <vector>

#include <seqan3/alphabet/nucleotide/nucleotide_base.hpp>
#include <seqan3/core/char_operations/transform.hpp>

// ------------------------------------------------------------------
// dna5
// ------------------------------------------------------------------

namespace seqan3
{

class rna5;

/*!\brief The five letter DNA alphabet of A,C,G,T and the unknown character N.
 * \ingroup nucleotide
 * \implements seqan3::nucleotide_alphabet
 * \implements seqan3::writable_alphabet
 * \if DEV \implements seqan3::detail::writable_constexpr_alphabet \endif
 * \implements seqan3::trivially_copyable
 * \implements seqan3::standard_layout
 * \implements std::regular
 *
 * \details
 * Note that you can assign 'U' as a character to dna5 and it will silently
 * be converted to 'T'.
 *
 * Like most alphabets, this alphabet cannot be initialised directly from its character representation.
 * Instead initialise/assign from the character literal or use the
 * function seqan3::dna5::assign_char().
 *
 *\include test/snippet/alphabet/nucleotide/dna5.cpp
 */
class dna5 : public nucleotide_base<dna5, 5>
{
private:
    //!\brief The base class.
    using base_t = nucleotide_base<dna5, 5>;

    //!\brief Befriend seqan3::nucleotide_base.
    friend base_t;
    //!\cond \brief Befriend seqan3::alphabet_base.
    friend base_t::base_t;
    //!\endcond
    //!\brief Befriend seqan3::rna5 so it can copy #char_to_rank.
    friend rna5;

public:
    /*!\name Constructors, destructor and assignment
     * \{
     */
    constexpr dna5()                          noexcept = default; //!< Defaulted.
    constexpr dna5(dna5 const &)              noexcept = default; //!< Defaulted.
    constexpr dna5(dna5 &&)                   noexcept = default; //!< Defaulted.
    constexpr dna5 & operator=(dna5 const &)  noexcept = default; //!< Defaulted.
    constexpr dna5 & operator=(dna5 &&)       noexcept = default; //!< Defaulted.
    ~dna5()                                   noexcept = default; //!< Defaulted.

    using base_t::base_t;

    //!\brief Allow implicit construction from dna/rna of the same size.
    template <std::same_as<rna5> t>    // Accept incomplete type
    constexpr dna5(t const & r) noexcept
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
        'N',
        'T'
    };

    //!\copydoc seqan3::dna4::char_to_rank
    static constexpr std::array<rank_type, 256> char_to_rank
    {
        [] () constexpr
        {
            std::array<rank_type, 256> ret{};

            // initialize with UNKNOWN (std::array::fill unfortunately not constexpr)
            for (auto & c : ret)
                c = 3; // == 'N'

            // reverse mapping for characters and their lowercase
            for (size_t rnk = 0u; rnk < alphabet_size; ++rnk)
            {
                ret[         rank_to_char[rnk] ] = rnk;
                ret[to_lower(rank_to_char[rnk])] = rnk;
            }

            // set U equal to T
            ret['U'] = ret['T']; ret['u'] = ret['t'];

            // iupac characters are implicitly "UNKNOWN"
            return ret;
        }()
    };

    //!\copydoc seqan3::dna4::complement_table
    static const std::array<dna5, alphabet_size> complement_table;
};

// ------------------------------------------------------------------
// containers
// ------------------------------------------------------------------

//!\brief Alias for an std::vector of seqan3::dna5.
//!\relates dna5
using dna5_vector = std::vector<dna5>;

// ------------------------------------------------------------------
// literals
// ------------------------------------------------------------------

/*!\name Literals
 * \{
 */

/*!\brief The seqan3::dna5 char literal.
 * \relates seqan3::dna5
 * \returns seqan3::dna5
 */
constexpr dna5 operator""_dna5(char const c) noexcept
{
    return dna5{}.assign_char(c);
}

/*!\brief The seqan3::dna5 string literal.
 * \relates seqan3::dna5
 * \returns seqan3::dna5_vector
 *
 * You can use this string literal to easily assign to dna5_vector:
 *
 * \include test/snippet/alphabet/nucleotide/dna5_literal.cpp
 *
 */
inline dna5_vector operator""_dna5(char const * s, std::size_t n)
{
    dna5_vector r;
    r.resize(n);

    for (size_t i = 0; i < n; ++i)
        r[i].assign_char(s[i]);

    return r;
}
//!\}

// ------------------------------------------------------------------
// dna5 (deferred definition)
// ------------------------------------------------------------------

constexpr std::array<dna5, dna5::alphabet_size> dna5::complement_table
{
    'T'_dna5,    // complement of 'A'_dna5
    'G'_dna5,    // complement of 'C'_dna5
    'C'_dna5,    // complement of 'G'_dna5
    'N'_dna5,    // complement of 'N'_dna5
    'A'_dna5     // complement of 'T'_dna5
};

} // namespace seqan3
