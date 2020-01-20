// -----------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2020, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2020, Knut Reinert & MPI für molekulare Genetik
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at: https://github.com/seqan/seqan3/blob/master/LICENSE.md
// -----------------------------------------------------------------------------------------------------

/*!\file
 * \author Sara Hetzel <sara.hetzel AT fu-berlin.de>
 * \brief Provides seqan3::aa27, container aliases and string literals.
 */

#pragma once

#include <vector>

#include <seqan3/alphabet/aminoacid/aminoacid_base.hpp>
#include <seqan3/alphabet/aminoacid/concept.hpp>
#include <seqan3/core/char_operations/transform.hpp>

namespace seqan3
{
/*!\brief The twenty-seven letter amino acid alphabet.
 * \ingroup aminoacid
 * \implements seqan3::aminoacid_alphabet
 * \implements seqan3::writable_alphabet
 * \if DEV \implements seqan3::detail::writable_constexpr_alphabet \endif
 * \implements seqan3::trivially_copyable
 * \implements seqan3::standard_layout
 * \implements std::regular
 *
 * \details
 * The alphabet consists of letters A, B, C, D, E, F, G, H, I, J, K, L, M, N, O, P, Q, R, S, T, U, V, W, X,
 * Y, Z, *
 *
 * Like most alphabets, this alphabet cannot be initialised directly from its character representation.
 * Instead initialise/assign from the character literal or use the
 * function seqan3::aa27::assign_char().
 *
 * \include test/snippet/alphabet/aminoacid/aa27_construction.cpp
 */

class aa27 : public aminoacid_base<aa27, 27>
{
private:
    //!\brief The base class.
    using base_t = aminoacid_base<aa27, 27>;

    //!\brief Befriend seqan3::nucleotide_base.
    friend base_t;
    //!\cond \brief Befriend seqan3::alphabet_base.
    friend base_t::base_t;
    //!\endcond

public:
    /*!\name Constructors, destructor and assignment
     * \{
     */
    constexpr aa27()                         noexcept = default; //!< Defaulted.
    constexpr aa27(aa27 const &)             noexcept = default; //!< Defaulted.
    constexpr aa27(aa27 &&)                  noexcept = default; //!< Defaulted.
    constexpr aa27 & operator=(aa27 const &) noexcept = default; //!< Defaulted.
    constexpr aa27 & operator=(aa27 &&)      noexcept = default; //!< Defaulted.
    ~aa27()                                  noexcept = default; //!< Defaulted.

    using base_t::base_t;
    //!\}

protected:
    //!\brief Value to char conversion table.
    static constexpr char_type rank_to_char[alphabet_size]
    {
        'A',
        'B',
        'C',
        'D',
        'E',
        'F',
        'G',
        'H',
        'I',
        'J',
        'K',
        'L',
        'M',
        'N',
        'O',
        'P',
        'Q',
        'R',
        'S',
        'T',
        'U',
        'V',
        'W',
        'X',
        'Y',
        'Z',
        '*'
    };

    //!\brief Char to value conversion table.
    static constexpr std::array<rank_type, 256> char_to_rank
    {
        [] () constexpr
        {
            std::array<rank_type, 256> ret{};

            // initialize with UNKNOWN (std::array::fill unfortunately not constexpr)
            for (auto & c : ret)
                c = 23; // value of 'X'

            // reverse mapping for characters and their lowercase
            for (rank_type rnk = 0u; rnk < alphabet_size; ++rnk)
            {
                ret[static_cast<rank_type>(         rank_to_char[rnk]) ] = rnk;
                ret[static_cast<rank_type>(to_lower(rank_to_char[rnk]))] = rnk;
            }

            return ret;
        }()
    };
};

} // namespace seqan3

// ------------------------------------------------------------------
// containers
// ------------------------------------------------------------------

namespace seqan3
{
//!\brief Alias for an std::vector of seqan3::aa27.
//!\relates aa27
using aa27_vector = std::vector<aa27>;

} // namespace seqan3

// ------------------------------------------------------------------
// literals
// ------------------------------------------------------------------

namespace seqan3
{

/*!\name Literals
 * \{
 */

/*!\brief The seqan3::aa27 char literal.
 * \param[in] c The character to assign.
 * \relates seqan3::aa27
 * \returns seqan3::aa27
 *
 * \include test/snippet/alphabet/aminoacid/aa27_char_literal.cpp
 *
 */
constexpr aa27 operator""_aa27(char const c) noexcept
{
    return aa27{}.assign_char(c);
}

/*!\brief The seqan3::aa27 string literal.
 * \param[in] s A pointer to the character string to assign.
 * \param[in] n The size of the character string to assign.
 * \relates seqan3::aa27
 * \returns seqan3::aa27_vector
 *
 * You can use this string literal to easily assign to aa27_vector:
 *
 * \include test/snippet/alphabet/aminoacid/aa27_literal.cpp
 *
 * \attention
 * All seqan3 literals are in the namespace seqan3!
 */

inline aa27_vector operator""_aa27(const char * s, std::size_t n)
{
    aa27_vector r;
    r.resize(n);

    for (size_t i = 0; i < n; ++i)
        r[i].assign_char(s[i]);

    return r;
}
//!\}

} // namespace seqan3
