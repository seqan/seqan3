// -----------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2020, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2020, Knut Reinert & MPI für molekulare Genetik
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at: https://github.com/seqan/seqan3/blob/master/LICENSE.md
// -----------------------------------------------------------------------------------------------------

/*!\file
 * \author Sara Hetzel <sara.hetzel AT fu-berlin.de>
 * \brief Provides seqan3::dna3bs, container aliases and string literals.
 */

#pragma once

#include <vector>

#include <seqan3/alphabet/nucleotide/nucleotide_base.hpp>
#include <seqan3/core/char_operations/transform.hpp>

// ------------------------------------------------------------------
// dna3bs
// ------------------------------------------------------------------

namespace seqan3
{
/*!\brief The three letter reduced DNA alphabet for bisulfite sequencing mode (A,G,T(=C)).
 * \ingroup nucleotide
 * \implements seqan3::nucleotide_alphabet
 * \implements seqan3::writable_alphabet
 * \if DEV \implements seqan3::detail::writable_constexpr_alphabet \endif
 * \implements seqan3::trivially_copyable
 * \implements seqan3::standard_layout
 * \implements std::regular
 *
 * \details
 * This alphabet represents a reduced version that can be used when dealing with bisulfite-converted data.
 * All 'C's are converted to a 'T' in order to allow comparison of normal sequences with
 * bisulfite-converted sequences.
 * For completeness, this nucleotide alphabet has a complement table, however, it is not recommended to use
 * it when dealing with bisulfite data because the complement of T is ambiguous in reads from bisulfite
 * sequencing. A 'T' can represent a true thymidine or an unmethylated 'C' that was converted into a 'T'.
 * Therefore, complementing a dna4bs sequence will further reduce the alphabet to only 'T' and 'A', thereby
 * loosing all information about 'G'. When working with bisulfite data, we recommend to create the reverse
 * complement of the dna4/5/15 range first and convert to dna3bs later. This avoids simplifying the data by automatically
 * setting 'A' as the complement of 'C'. As an example: The sequence 'ACGTGC' in dna4 would be 'ATGTGT' in dna3bs.
 * The complement of this dna3bs sequence would be 'TATATA', however when complementing the dna4 sequence first
 * and afterwards transforming it into dna3bs, it would be 'TGTATG' which preserves more information from the original sequence.
 *
 * Like most alphabets, this alphabet cannot be initialised directly from its character representation.
 * Instead initialise/assign from the character literal or use the
 * function seqan3::dna3bs::assign_char().
 *
 * \include test/snippet/alphabet/nucleotide/dna3bs.cpp
 */
class dna3bs : public nucleotide_base<dna3bs, 3>
{
private:
    //!\brief The base class.
    using base_t = nucleotide_base<dna3bs, 3>;

    //!\brief Befriend seqan3::nucleotide_base.
    friend base_t;
    //!\cond \brief Befriend seqan3::alphabet_base.
    friend base_t::base_t;
    //!\endcond

public:
    /*!\name Constructors, destructor and assignment
     * \{
     */
    constexpr dna3bs()                            noexcept = default; //!< Defaulted.
    constexpr dna3bs(dna3bs const &)              noexcept = default; //!< Defaulted.
    constexpr dna3bs(dna3bs &&)                   noexcept = default; //!< Defaulted.
    constexpr dna3bs & operator=(dna3bs const &)  noexcept = default; //!< Defaulted.
    constexpr dna3bs & operator=(dna3bs &&)       noexcept = default; //!< Defaulted.
    ~dna3bs()                                     noexcept = default; //!< Defaulted.

    using base_t::base_t;
    //!\}

protected:
    //!\privatesection

    //!\brief Value to char conversion table.
    static constexpr char_type rank_to_char[alphabet_size]
    {
        'A',
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
            for (size_t rnk = 0u; rnk < alphabet_size; ++rnk)
            {
                ret[         rank_to_char[rnk] ] = rnk;
                ret[to_lower(rank_to_char[rnk])] = rnk;
            }

            // set C and U equal to T
            ret['C'] = ret['T']; ret['c'] = ret['t'];
            ret['U'] = ret['T']; ret['u'] = ret['t'];

            // iupac characters get special treatment, because there is no N
            ret['R'] = ret['A']; ret['r'] = ret['A']; // A or G becomes A
            ret['Y'] = ret['T']; ret['y'] = ret['T']; // C or T becomes T
            ret['S'] = ret['T']; ret['s'] = ret['T']; // C or G becomes T
            ret['W'] = ret['A']; ret['w'] = ret['A']; // A or T becomes A
            ret['K'] = ret['G']; ret['k'] = ret['G']; // G or T becomes G
            ret['M'] = ret['A']; ret['m'] = ret['A']; // A or C becomes A
            ret['B'] = ret['T']; ret['b'] = ret['T']; // C or G or T becomes T
            ret['D'] = ret['A']; ret['d'] = ret['A']; // A or G or T becomes A
            ret['H'] = ret['A']; ret['h'] = ret['A']; // A or C or T becomes A
            ret['V'] = ret['A']; ret['v'] = ret['A']; // A or C or G  becomes A

            return ret;
        }()
    };

    //!\brief The complement table.
    static const std::array<dna3bs, alphabet_size> complement_table;
};

// ------------------------------------------------------------------
// containers
// ------------------------------------------------------------------

//!\brief Alias for an std::vector of seqan3::dna3bs.
//!\relates dna3bs
using dna3bs_vector = std::vector<dna3bs>;

// ------------------------------------------------------------------
// literals
// ------------------------------------------------------------------

/*!\name Literals
 * \{
 */

/*!\brief The seqan3::dna3bs char literal.
 * \relates seqan3::dna3bs
 * \returns seqan3::dna3bs
 */
constexpr dna3bs operator""_dna3bs(char const c) noexcept
{
    return dna3bs{}.assign_char(c);
}

/*!\brief The seqan3::dna3bs string literal.
 * \relates seqan3::dna3bs
 * \returns seqan3::dna3bs_vector
 *
 * You can use this string literal to easily assign to dna3bs_vector:
 *
 * \include test/snippet/alphabet/nucleotide/dna3bs_literal.cpp
 *
 */
inline dna3bs_vector operator""_dna3bs(char const * s, std::size_t n)
{
    dna3bs_vector r;
    r.resize(n);

    for (size_t i = 0; i < n; ++i)
        r[i].assign_char(s[i]);

    return r;
}
//!\}

// ------------------------------------------------------------------
// dna3bs (deferred definition)
// ------------------------------------------------------------------

constexpr std::array<dna3bs, dna3bs::alphabet_size> dna3bs::complement_table
{
    'T'_dna3bs,    // complement of 'A'_dna3bs
    'T'_dna3bs,    // complement of 'G'_dna3bs
    'A'_dna3bs     // complement of 'T'_dna3bs
};

} // namespace seqan3
