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
 * \brief Contains seqan3::phred42 quality scores.
 */

#pragma once

#include <seqan3/alphabet/quality/quality_base.hpp>

// ------------------------------------------------------------------
// phred42
// ------------------------------------------------------------------

namespace seqan3
{

/*!\brief Quality type for traditional Sanger and modern Illumina Phred scores (typical range).
 * \implements seqan3::quality_concept
 * \implements seqan3::detail::constexpr_alphabet_concept
 * \implements seqan3::trivially_copyable_concept
 * \implements seqan3::standard_layout_concept
 *
 * \ingroup quality
 *
 * \details
 *
 * The phred42 quality alphabet represents the zero-based phred score range
 * [0..41] mapped to the consecutive ASCII range ['!' .. 'J']. It therefore can
 * represent the Illumina 1.8+ standard and the original Sanger score. If you
 * intend to use phred scores exceeding 41, use the larger score type, namely
 * seqan3::phred63, otherwise on construction exceeding scores are mapped to 41.
 *
 * \snippet test/snippet/alphabet/quality/phred42.cpp general
 */
class phred42 : public quality_base<phred42, 42>
{
private:
    //!\brief The base class.
    using base_t = quality_base<phred42, 42>;

    //!\brief Befriend seqan3::quality_base.
    friend base_t;
    //!\cond \brief Befriend seqan3::alphabet_base.
    friend base_t::base_t;
    //!\endcond

public:
    /*!\name Constructors, destructor and assignment
     * \{
     */
    constexpr phred42() noexcept : base_t{} {}
    constexpr phred42(phred42 const &) = default;
    constexpr phred42(phred42 &&) = default;
    constexpr phred42 & operator=(phred42 const &) = default;
    constexpr phred42 & operator=(phred42 &&) = default;
    ~phred42() = default;

    //!\brief Construct from phred value.
    constexpr phred42(phred_type const p) : base_t{p} {}

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
/*!\brief The seqan3::phred42 char literal.
 * \relates seqan3::phred42
 * \returns seqan3::phred42
 */
constexpr phred42 operator""_phred42(char const c) noexcept
{
    return phred42{}.assign_char(c);
}

/*!\brief The seqan3::phred42 string literal.
 * \param[in] s A pointer to the character sequence to assign from.
 * \param[in] n The length of the character sequence to assign from.
 * \relates seqan3::phred42
 * \returns seqan3::std::vector<seqan3::phred42>
 *
 * You can use this string literal to easily assign to std::vector<seqan3::phred42>:
 *
 * \include test/snippet/alphabet/quality/phred42_literal.cpp
 */
inline std::vector<phred42> operator""_phred42(char const * s, std::size_t n)
{
    std::vector<phred42> r;
    r.resize(n);

    for (size_t i = 0; i < n; ++i)
        r[i].assign_char(s[i]);

    return r;
}
//!\}

} // namespace seqan3
