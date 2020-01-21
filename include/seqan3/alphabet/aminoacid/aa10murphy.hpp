// -----------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2020, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2020, Knut Reinert & MPI für molekulare Genetik
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at: https://github.com/seqan/seqan3/blob/master/LICENSE.md
// -----------------------------------------------------------------------------------------------------

/*!\file
 * \author Sara Hetzel <sara.hetzel AT fu-berlin.de>
 * \brief Provides seqan3::aa10murphy, container aliases and string literals.
 */

#pragma once

#include <vector>

#include <seqan3/alphabet/aminoacid/aminoacid_base.hpp>
#include <seqan3/alphabet/aminoacid/concept.hpp>
#include <seqan3/core/char_operations/transform.hpp>

namespace seqan3
{

/*!\brief The reduced Murphy amino acid alphabet.
 * \ingroup aminoacid
 * \implements seqan3::aminoacid_alphabet
 * \implements seqan3::writable_alphabet
 * \if DEV \implements seqan3::detail::writable_constexpr_alphabet \endif
 * \implements seqan3::trivially_copyable
 * \implements seqan3::standard_layout
 * \implements std::regular
 *
 * \details
 * The alphabet consists of letters A, B, C, F, G, H, I, K, P, S
 * B represents charged/polar residues (E,D,N,Q).
 * C represents cystein and the species-specific amino acid Selenocysteine.
 * F represents amino acids with large and mainly hydrophobic aromatic side chains (F,W,Y).
 * I represents large hydrophobes (L,V,I,M).
 * K represents long-chain positively charged residues (K,R) and the species-specific amino acid Pyrrolysine.
 * S represents alcohols (S,T) and unknown.
 * A, G, H and P do not represent any other amino acids other than themselves.
 *
 * This alphabet allows to reduce the aminoacid alphabet size to 10 but is still able to recognize and represent
 * folding of all proteins. Amino acids are grouped together based on similar physical and chemical properties.
 *
 * *Note:* Letters which belong in the extended alphabet will be automatically converted.
 * Terminator characters are converted to F, because the most commonly occurring stop codon in higher eukaryotes
 * is UGA <sup>2</sup>. This is most similar to a Tryptophan which in this alphabet
 * gets converted to Phenylalanine. Anything unknown is converted to S.
 *
 * |Input Letter  |Converts to    |
 * |--------------|---------------|
 * |D             |B<sup>1</sup>  |
 * |E             |B<sup>1</sup>  |
 * |J             |I<sup>1</sup>  |
 * |L             |I<sup>1</sup>  |
 * |M             |I<sup>1</sup>  |
 * |N             |B<sup>1</sup>  |
 * |O             |K<sup>1</sup>  |
 * |Q             |B<sup>1</sup>  |
 * |R             |K<sup>1</sup>  |
 * |T             |S<sup>1</sup>  |
 * |U             |C<sup>1</sup>  |
 * |V             |I<sup>1</sup>  |
 * |W             |F<sup>1</sup>  |
 * |Y             |F<sup>1</sup>  |
 * |Z             |B<sup>1</sup>  |
 * |X (Unknown)   |S<sup>1</sup>  |
 * |* (Terminator)|F<sup>1,2</sup>|
 *
 * <sup><b>1</b></sup>L. R. Murphy, A. Wallqvist, and R. M. Levy. Simplified amino acid alphabets for protein
 * fold recognition and implications for folding. Protein Eng., 13(3):149–152, Mar 2000.\n
 * <sup><b>2</b></sup>Trotta, E. (2016). Selective forces and mutational biases drive stop codon usage
 * in the human genome: a comparison with sense codon usage.
 * BMC Genomics, 17, 366. https://doi.org/10.1186/s12864-016-2692-4
 *
 * \include test/snippet/alphabet/aminoacid/aa10murphy.cpp
 */
class aa10murphy : public aminoacid_base<aa10murphy, 10>
{
private:
    //!\brief The base class.
    using base_t = aminoacid_base<aa10murphy, 10>;

    //!\brief Befriend seqan3::aminoacid_base.
    friend base_t;
    //!\cond \brief Befriend seqan3::alphabet_base.
    friend base_t::base_t;
    //!\endcond

public:
    /*!\name Constructors, destructor and assignment
     * \{
     */
    constexpr aa10murphy()                               noexcept = default; //!< Defaulted.
    constexpr aa10murphy(aa10murphy const &)             noexcept = default; //!< Defaulted.
    constexpr aa10murphy(aa10murphy &&)                  noexcept = default; //!< Defaulted.
    constexpr aa10murphy & operator=(aa10murphy const &) noexcept = default; //!< Defaulted.
    constexpr aa10murphy & operator=(aa10murphy &&)      noexcept = default; //!< Defaulted.
    ~aa10murphy()                                        noexcept = default; //!< Defaulted.

    //!\brief Inherit the base class's Constructors.
    using base_t::base_t;
    //!\}

protected:
    //!\brief Value to char conversion table.
    static constexpr char_type rank_to_char[alphabet_size]
    {
        'A',
        'B',
        'C',
        'F',
        'G',
        'H',
        'I',
        'K',
        'P',
        'S',
    };

    //!\brief Char to value conversion table.
    static constexpr std::array<rank_type, 256> char_to_rank
    {
        [] () constexpr
        {
            std::array<rank_type, 256> ret{};

            // initialize with UNKNOWN (std::array::fill unfortunately not constexpr)
            for (auto & c : ret)
                c = 9; // value of 'S', because that appears most frequently

            // reverse mapping for characters and their lowercase
            for (rank_type rnk = 0u; rnk < alphabet_size; ++rnk)
            {
                ret[static_cast<rank_type>(         rank_to_char[rnk]) ] = rnk;
                ret[static_cast<rank_type>(to_lower(rank_to_char[rnk]))] = rnk;
            }

            ret['D'] = ret['B']; ret['d'] = ret['B']; // Convert D to B (either D/N).
            ret['E'] = ret['B']; ret['e'] = ret['B']; // Convert E to B (either D/N).
            ret['J'] = ret['I']; ret['j'] = ret['I']; // Convert J (either I/L) to I.
            ret['L'] = ret['I']; ret['l'] = ret['I']; // Convert L to I.
            ret['M'] = ret['I']; ret['m'] = ret['I']; // Convert M to I.
            ret['N'] = ret['B']; ret['n'] = ret['B']; // Convert N to B (either D/N).
            ret['O'] = ret['K']; ret['o'] = ret['K']; // Convert Pyrrolysine to K.
            ret['Q'] = ret['B']; ret['q'] = ret['B']; // Convert Q to B (either D/N).
            ret['R'] = ret['K']; ret['r'] = ret['K']; // Convert R to K.
            ret['T'] = ret['S']; ret['t'] = ret['S']; // Convert T to S.
            ret['U'] = ret['C']; ret['u'] = ret['C']; // Convert Selenocysteine to C.
            ret['V'] = ret['I']; ret['v'] = ret['I']; // Convert V to I.
            ret['W'] = ret['F']; ret['w'] = ret['F']; // Convert W to F.
            ret['X'] = ret['S']; ret['x'] = ret['S']; // Convert unknown amino acids to Serine.
            ret['Y'] = ret['F']; ret['y'] = ret['F']; // Convert Y to F.
            ret['Z'] = ret['B']; ret['z'] = ret['B']; // Convert Z (either E/Q) to B (either D/N).
            ret['*'] = ret['F']; // The most common stop codon is UGA. This is most similar to a Tryptophan which in this alphabet gets converted to Phenylalanine.
            return ret;
        }()
    };
};

// ------------------------------------------------------------------
// containers
// ------------------------------------------------------------------

//!\brief Alias for an std::vector of seqan3::aa10murphy.
//!\relates aa10murphy
using aa10murphy_vector = std::vector<aa10murphy>;

// ------------------------------------------------------------------
// literals
// ------------------------------------------------------------------

/*!\name Literals
 * \{
 */

/*!\brief The seqan3::aa10murphy char literal.
 * \param[in] c The character to assign.
 * \relates seqan3::aa10murphy
 * \returns seqan3::aa10murphy
 */
constexpr aa10murphy operator""_aa10murphy(char const c) noexcept
{
    return aa10murphy{}.assign_char(c);
}

/*!\brief The seqan3::aa10murphy  string literal.
 * \param[in] s A pointer to the character string to assign.
 * \param[in] n The size of the character string to assign.
 * \relates seqan3::aa10murphy
 * \returns seqan3::aa10murphy_vector
 *
 * You can use this string literal to easily assign to aa10murphy_vector:
 *
 * \attention
 * All seqan3 literals are in the namespace seqan3!
 */

inline aa10murphy_vector operator""_aa10murphy(const char * s, std::size_t n)
{
    aa10murphy_vector r;
    r.resize(n);

    for (size_t i = 0; i < n; ++i)
        r[i].assign_char(s[i]);

    return r;
}
//!\}

} // namespace seqan3
