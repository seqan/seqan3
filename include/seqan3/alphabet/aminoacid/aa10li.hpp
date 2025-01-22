// SPDX-FileCopyrightText: 2006-2025 Knut Reinert & Freie Universität Berlin
// SPDX-FileCopyrightText: 2016-2025 Knut Reinert & MPI für molekulare Genetik
// SPDX-License-Identifier: BSD-3-Clause

/*!\file
 * \author Sara Hetzel <sara.hetzel AT fu-berlin.de>
 * \brief Provides seqan3::aa10li, container aliases and string literals.
 */

#pragma once

#include <vector>

#include <seqan3/alphabet/aminoacid/aminoacid_base.hpp>
#include <seqan3/alphabet/aminoacid/concept.hpp>
#include <seqan3/utility/char_operations/transform.hpp>

namespace seqan3
{
/*!\brief The reduced Li amino acid alphabet.
 * \ingroup alphabet_aminoacid
 * \implements seqan3::aminoacid_alphabet
 * \implements seqan3::writable_alphabet
 * \implements seqan3::detail::writable_constexpr_alphabet
 * \implements seqan3::trivially_copyable
 * \implements seqan3::standard_layout
 * \implements std::regular
 *
 * \details
 * The alphabet consists of letters A, B, C, F, G, H, I, J, K, P
 * A represents hydrophilic and alocohol residues (A,S,T).
 * B represents charged/polar residues (B,D,E,Q,Z).
 * C represents cystein and the species-specific amino acid Selenocysteine.
 * F represents amino acids with aromatic residues (F,W,Y).
 * H represents a group of hydrophobic residues (H,N).
 * I represents a group of large hydrophobic residues (I,V).
 * J represents a group of large hydrophobic residues (J,L,M).
 * K represents long-chain positively charged residues (K,R) and the species-specific amino acid Pyrrolysine.
 * G and P do not represent any other amino acids other than themselves.
 *
 * This alphabet allows to reduce the aminoacid alphabet size to 10 but is still able to recognize and represent
 * folding of all proteins. Amino acids are grouped together based on residues.
 *
 * *Note:* Letters which belong in the extended alphabet will be automatically converted.
 * Terminator characters are converted to F, because the most commonly occurring stop codon in higher eukaryotes
 * is UGA <sup>2</sup>. This is most similar to a Tryptophan which in this alphabet
 * gets converted to Phenylalanine. Anything unknown is converted to A.
 *
 * |Input Letter  |Converts to    |
 * |--------------|---------------|
 * |D             |B<sup>1</sup>  |
 * |E             |B<sup>1</sup>  |
 * |L             |J<sup>1</sup>  |
 * |M             |J<sup>1</sup>  |
 * |N             |H<sup>1</sup>  |
 * |O             |K<sup>1</sup>  |
 * |Q             |B<sup>1</sup>  |
 * |R             |K<sup>1</sup>  |
 * |S             |A<sup>1</sup>  |
 * |T             |A<sup>1</sup>  |
 * |U             |C<sup>1</sup>  |
 * |V             |I<sup>1</sup>  |
 * |W             |F<sup>1</sup>  |
 * |Y             |F<sup>1</sup>  |
 * |Z             |B<sup>1</sup>  |
 * |X (Unknown)   |A<sup>1</sup>  |
 * |* (Terminator)|F<sup>1,2</sup>|
 *
 * <sup><b>1</b></sup>T. Li, K. Fan, J. Wang, and W. Wang. Reduction of protein sequence complexity by
 * residue grouping. Protein Eng., 16(5):323–330, May 2003.\n
 * <sup><b>2</b></sup>Trotta, E. (2016). Selective forces and mutational biases drive stop codon usage
 * in the human genome: a comparison with sense codon usage.
 * BMC Genomics, 17, 366. https://doi.org/10.1186/s12864-016-2692-4
 *
 * \include test/snippet/alphabet/aminoacid/aa10li.cpp
 *
 * \stableapi{Since version 3.1.}
 */
class aa10li : public aminoacid_base<aa10li, 10>
{
private:
    //!\brief The base class.
    using base_t = aminoacid_base<aa10li, 10>;

    //!\brief Befriend seqan3::aminoacid_base.
    friend base_t;
    //!\cond
    //!\brief Befriend seqan3::alphabet_base.
    friend base_t::base_t;
    //!\endcond

public:
    /*!\name Constructors, destructor and assignment
     * \{
     */
    constexpr aa10li() noexcept = default;                           //!< Defaulted.
    constexpr aa10li(aa10li const &) noexcept = default;             //!< Defaulted.
    constexpr aa10li(aa10li &&) noexcept = default;                  //!< Defaulted.
    constexpr aa10li & operator=(aa10li const &) noexcept = default; //!< Defaulted.
    constexpr aa10li & operator=(aa10li &&) noexcept = default;      //!< Defaulted.
    ~aa10li() noexcept = default;                                    //!< Defaulted.

    //!\brief Inherit the base class's Constructors.
    using base_t::base_t;
    //!\}

private:
    //!\copydoc seqan3::aa27::rank_to_char_table
    static constexpr char_type rank_to_char_table[alphabet_size]{'A', 'B', 'C', 'F', 'G', 'H', 'I', 'J', 'K', 'P'};

    //!\copydoc seqan3::aa27::char_to_rank
    static constexpr rank_type char_to_rank(char_type const chr)
    {
        using index_t = std::make_unsigned_t<char_type>;
        return char_to_rank_table[static_cast<index_t>(chr)];
    }

    //!\copydoc seqan3::aa27::rank_to_char
    static constexpr char_type rank_to_char(rank_type const rank)
    {
        return rank_to_char_table[rank];
    }

    //!\copydoc seqan3::aa27::char_to_rank_table
    static constexpr std::array<rank_type, 256> char_to_rank_table{
        []() constexpr
        {
            std::array<rank_type, 256> ret{};

            // initialize with 'A' because S appears most frequently and gets converted to A in this alphabet
            ret.fill(0u); // Value-initialisation of std::array does usually initialise. `fill` is explicit.

            // reverse mapping for characters and their lowercase
            for (rank_type rnk = 0u; rnk < alphabet_size; ++rnk)
            {
                ret[static_cast<rank_type>(rank_to_char_table[rnk])] = rnk;
                ret[static_cast<rank_type>(to_lower(rank_to_char_table[rnk]))] = rnk;
            }

            ret['D'] = ret['B'];
            ret['d'] = ret['B']; // Convert D to B (either D/N).
            ret['E'] = ret['B'];
            ret['e'] = ret['B']; // Convert E to B (either D/N).
            ret['L'] = ret['J'];
            ret['l'] = ret['J']; // Convert L to J (either I/L).
            ret['M'] = ret['J'];
            ret['m'] = ret['J']; // Convert M to J (either I/L).
            ret['N'] = ret['H'];
            ret['n'] = ret['H']; // Convert N to H.
            ret['O'] = ret['K'];
            ret['o'] = ret['K']; // Convert Pyrrolysine to K.
            ret['Q'] = ret['B'];
            ret['q'] = ret['B']; // Convert Q to B (either D/N).
            ret['R'] = ret['K'];
            ret['r'] = ret['K']; // Convert R to K.
            ret['S'] = ret['A'];
            ret['s'] = ret['A']; // Convert S to A.
            ret['T'] = ret['A'];
            ret['t'] = ret['A']; // Convert T to A.
            ret['U'] = ret['C'];
            ret['u'] = ret['C']; // Convert Selenocysteine to C.
            ret['V'] = ret['I'];
            ret['v'] = ret['I']; // Convert V to I.
            ret['W'] = ret['F'];
            ret['w'] = ret['F']; // Convert W to F.
            ret['X'] = ret['A'];
            ret['x'] = ret['A']; // Convert unknown amino acids to Alanine.
            ret['Y'] = ret['F'];
            ret['y'] = ret['F']; // Convert Y to F.
            ret['Z'] = ret['B'];
            ret['z'] = ret['B']; // Convert Z (either E/Q) to B (either D/N).
            ret['*'] = ret['F']; // The most common stop codon is UGA. This is most similar to a Tryptophan which in
                                 // this alphabet gets converted to Phenylalanine.

            return ret;
        }()};
};

// ------------------------------------------------------------------
// containers
// ------------------------------------------------------------------

/*!\brief Alias for a std::vector of seqan3::aa10li.
 * \relates aa10li
 *
 * \stableapi{Since version 3.1.}
 */
using aa10li_vector = std::vector<aa10li>;

// ------------------------------------------------------------------
// literals
// ------------------------------------------------------------------
inline namespace literals
{

/*!\name Literals
 * \{
 */
/*!\brief The seqan3::aa10li char literal.
 * \param[in] c The character to assign.
 * \relates seqan3::aa10li
 * \returns seqan3::aa10li
 *
 * You can use this char literal to assign a seqan3::aa10li character:
 * \include test/snippet/alphabet/aminoacid/aa10li_char_literal.cpp
 *
 * \stableapi{Since version 3.1.}
 */
constexpr aa10li operator""_aa10li(char const c) noexcept
{
    return aa10li{}.assign_char(c);
}

/*!\brief The seqan3::aa10li string literal.
 * \param[in] s A pointer to the character string to assign.
 * \param[in] n The size of the character string to assign.
 * \relates seqan3::aa10li
 * \returns seqan3::aa10li_vector
 *
 * You can use this string literal to easily assign to aa10li_vector:
 * \include test/snippet/alphabet/aminoacid/aa10li_literal.cpp
 *
 * \stableapi{Since version 3.1.}
 */
SEQAN3_WORKAROUND_LITERAL aa10li_vector operator""_aa10li(char const * const s, size_t const n)
{
    aa10li_vector r;
    r.resize(n);

    for (size_t i = 0; i < n; ++i)
        r[i].assign_char(s[i]);

    return r;
}
//!\}

} // namespace literals

} // namespace seqan3
