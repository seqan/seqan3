// ============================================================================
//                 SeqAn - The Library for Sequence Analysis
// ============================================================================
//
// Copyright (c) 2006-2018, Knut Reinert & Freie Universitaet Berlin
// Copyright (c) 2016-2018, Knut Reinert & MPI Molekulare Genetik
// All rights reserved.
//
// Redistribution and use in source and binary forms, with or without
// modification, are permitted provided that the following conditions are met:
//
//     * Redistributions of source code must retain the above copyright
//       notice, this list of conditions and the following disclaimer.
//     * Redistributions in binary form must reproduce the above copyright
//       notice, this list of conditions and the following disclaimer in the
//       documentation and/or other materials provided with the distribution.
//     * Neither the name of Knut Reinert or the FU Berlin nor the names of
//       its contributors may be used to endorse or promote products derived
//       from this software without specific prior written permission.
//
// THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
// AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
// IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
// ARE DISCLAIMED. IN NO EVENT SHALL KNUT REINERT OR THE FU BERLIN BE LIABLE
// FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL
// DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR
// SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER
// CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT
// LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY
// OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH
// DAMAGE.
//
// ============================================================================

/*!\file
 * \author Marie Hoffmann <marie.hoffmann AT fu-berlin.de>
 * \brief Contains seqan3::phred63 quality scores.
 */

#pragma once

#include <seqan3/alphabet/quality/quality_base.hpp>

// ------------------------------------------------------------------
// phred63
// ------------------------------------------------------------------

namespace seqan3
{

/*!\brief Quality type for traditional Sanger and modern Illumina Phred scores (full range).
 * \implements seqan3::quality_concept
 * \implements seqan3::detail::constexpr_alphabet_concept
 * \implements seqan3::trivially_copyable_concept
 * \implements seqan3::standard_layout_concept
 *
 * \ingroup quality
 *
 * \details
 *
 * The phred63 quality alphabet represents the zero-based phred score range
 * [0..62] mapped to the ASCII range ['!' .. '~']. It represents the Sanger and
 * Illumina 1.8+ standard beyond the typical range of raw reads (0 to 41).
 *
 * \snippet test/snippet/alphabet/quality/phred63.cpp general
 */
class phred63 : public quality_base<phred63, 63>
{
private:
    //!\brief The base class.
    using base_t = quality_base<phred63, 63>;

    //!\brief Befriend seqan3::quality_base.
    friend base_t;
    //!\cond \brief Befriend seqan3::alphabet_base.
    friend base_t::base_t;
    //!\endcond

public:
    /*!\name Constructors, destructor and assignment
     * \{
     */
    constexpr phred63() noexcept : base_t{} {}
    constexpr phred63(phred63 const &) = default;
    constexpr phred63(phred63 &&) = default;
    constexpr phred63 & operator=(phred63 const &) = default;
    constexpr phred63 & operator=(phred63 &&) = default;
    ~phred63() = default;

    //!\brief Construct from phred value.
    constexpr phred63(phred_type const p) : base_t{p} {}

    // Inherit converting constructor
    using base_t::base_t;
    //!\}

    /*!\name Member variables.
     * \{
     */
    //!\brief The projection offset between phred and rank score representation.
    static constexpr phred_type offset_phred{0};

    //!\brief The projection offset between char and rank score representation.
    static constexpr char_type offset_char{'!'};
    //!\}
};

/*!\name Literals
 * \{
 */
/*!\brief The seqan3::phred63 char literal.
 * \relates seqan3::phred63
 * \returns seqan3::phred63
 */
constexpr phred63 operator""_phred63(char const c) noexcept
{
    return phred63{}.assign_char(c);
}

/*!\brief The seqan3::phred63 string literal.
 * \param[in] s A pointer to the character sequence to assign from.
 * \param[in] n The length of the character sequence to assign from.
 * \relates seqan3::phred63
 * \returns seqan3::std::vector<seqan3::phred63>
 *
 * You can use this string literal to easily assign to std::vector<seqan3::phred63>:
 *
 * \include test/snippet/alphabet/quality/phred63_literal.cpp
 */
inline std::vector<phred63> operator""_phred63(char const * s, std::size_t n)
{
    std::vector<phred63> r;
    r.resize(n);

    for (size_t i = 0; i < n; ++i)
        r[i].assign_char(s[i]);

    return r;
}
//!\}

} // namespace seqan3
