// -----------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2019, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2019, Knut Reinert & MPI für molekulare Genetik
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at: https://github.com/seqan/seqan3/blob/master/LICENSE.md
// -----------------------------------------------------------------------------------------------------

/*!\file
 * \author Hannes Hauswedell <hannes.hauswedell AT fu-berlin.de>
 * \brief Contains seqan3::dna4, container aliases and string literals.
 */

#pragma once

#include <vector>

#include <seqan3/alphabet/nucleotide/nucleotide_base.hpp>
#include <seqan3/io/stream/char_operations.hpp>

// ------------------------------------------------------------------
// dna4
// ------------------------------------------------------------------

namespace seqan3
{

class rna4;

/*!\brief The four letter DNA alphabet of A,C,G,T.
 * \ingroup nucleotide
 * \implements seqan3::NucleotideAlphabet
 * \implements seqan3::WritableAlphabet
 * \if DEV \implements seqan3::detail::WritableConstexprAlphabet \endif
 * \implements seqan3::TriviallyCopyable
 * \implements seqan3::StandardLayout
 * \implements std::Regular
 *
 * \details
 * Note that you can assign 'U' as a character to dna4 and it will silently
 * be converted to 'T'.
 *
 * The alphabet may be brace initialized from the static letter members. Note that you cannot
 * assign the alphabet by using letters of type `char`, but you instead have to use the
 * function seqan3::dna4::assign_char().
 *
 * \snippet test/snippet/alphabet/nucleotide/dna4.cpp code
 */
class dna4 : public nucleotide_base<dna4, 4>
{
private:
    //!\brief The base class.
    using base_t = nucleotide_base<dna4, 4>;

    //!\brief Befriend seqan3::nucleotide_base.
    friend base_t;
    //!\cond \brief Befriend seqan3::alphabet_base.
    friend base_t::base_t;
    //!\endcond
    //!\brief Befriend seqan3::rna4 so it can copy #char_to_rank.
    friend rna4;

public:
    /*!\name Constructors, destructor and assignment
     * \{
     */
    constexpr dna4()                          noexcept = default; //!< Defaulted.
    constexpr dna4(dna4 const &)              noexcept = default; //!< Defaulted.
    constexpr dna4(dna4 &&)                   noexcept = default; //!< Defaulted.
    constexpr dna4 & operator=(dna4 const &)  noexcept = default; //!< Defaulted.
    constexpr dna4 & operator=(dna4 &&)       noexcept = default; //!< Defaulted.
    ~dna4()                                   noexcept = default; //!< Defaulted.

    using base_t::base_t;

    //!\brief Allow implicit construction from dna/rna of the same size.
    template <std::Same<rna4> t>    // Accept incomplete type
    constexpr dna4(t const & r) noexcept
    {
        assign_rank(r.to_rank());
    }
    //!\}

protected:
    //!\privatesection

    //!\brief Value to char conversion table.
    static constexpr char_type rank_to_char[value_size]
    {
        'A',
        'C',
        'G',
        'T'
    };

    //!\brief Char to value conversion table.
    static constexpr std::array<rank_type, 256> char_to_rank
    {
        [] () constexpr
        {
            std::array<rank_type, 256> ret{};

            // reverse mapping for characters and their lowercase
            for (size_t rnk = 0u; rnk < value_size; ++rnk)
            {
                ret[         rank_to_char[rnk] ] = rnk;
                ret[to_lower(rank_to_char[rnk])] = rnk;
            }

            // set U equal to T
            ret['U'] = ret['T']; ret['u'] = ret['t'];

            // iupac characters get special treatment, because there is no N
            ret['R'] = ret['A']; ret['r'] = ret['A']; // or G
            ret['Y'] = ret['C']; ret['y'] = ret['C']; // or T
            ret['S'] = ret['C']; ret['s'] = ret['C']; // or G
            ret['W'] = ret['A']; ret['w'] = ret['A']; // or T
            ret['K'] = ret['G']; ret['k'] = ret['G']; // or T
            ret['M'] = ret['A']; ret['m'] = ret['A']; // or T
            ret['B'] = ret['C']; ret['b'] = ret['C']; // or G or T
            ret['D'] = ret['A']; ret['d'] = ret['A']; // or G or T
            ret['H'] = ret['A']; ret['h'] = ret['A']; // or C or T
            ret['V'] = ret['A']; ret['v'] = ret['A']; // or C or G

            return ret;
        }()
    };

    //!\brief The complement table.
    static const std::array<dna4, value_size> complement_table;
};

// ------------------------------------------------------------------
// containers
// ------------------------------------------------------------------

//!\brief Alias for an std::vector of seqan3::dna4.
//!\relates dna4
using dna4_vector = std::vector<dna4>;

// ------------------------------------------------------------------
// literals
// ------------------------------------------------------------------

/*!\name Literals
 * \{
 */

/*!\brief The seqan3::dna4 char literal.
 * \relates seqan3::dna4
 * \returns seqan3::dna4
 */
constexpr dna4 operator""_dna4(char const c) noexcept
{
    return dna4{}.assign_char(c);
}

/*!\brief The seqan3::dna4 string literal.
 * \relates seqan3::dna4
 * \returns seqan3::dna4_vector
 *
 * You can use this string literal to easily assign to dna4_vector:
 *
 * \snippet test/snippet/alphabet/nucleotide/dna4.cpp operator""_dna4
 *
 */
inline dna4_vector operator""_dna4(char const * s, std::size_t n)
{
    dna4_vector r;
    r.resize(n);

    for (size_t i = 0; i < n; ++i)
        r[i].assign_char(s[i]);

    return r;
}
//!\}

// ------------------------------------------------------------------
// dna4 (deferred definition)
// ------------------------------------------------------------------

constexpr std::array<dna4, dna4::value_size> dna4::complement_table
{
    'T'_dna4,    // complement of 'A'_dna4
    'G'_dna4,    // complement of 'C'_dna4
    'C'_dna4,    // complement of 'G'_dna4
    'A'_dna4     // complement of 'T'_dna4
};

} // namespace seqan3
