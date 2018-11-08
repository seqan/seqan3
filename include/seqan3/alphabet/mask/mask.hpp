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
 * \author Joshua Kim <joshua.kim AT fu-berlin.de>
 * \brief Create a mask composition which can be applied with another alphabet.
 */

#pragma once

#include <cassert>
#include <seqan3/alphabet/concept_pre.hpp>
#include <seqan3/alphabet/detail/alphabet_base.hpp>

namespace seqan3
{
/*!\brief Implementation of a masked alphabet to be used for cartesian compositions.
 * \ingroup mask
 * \implements seqan3::semi_alphabet_concept
 * \implements seqan3::detail::semi_constexpr_alphabet_concept
 * \implements seqan3::trivially_copyable_concept
 * \implements seqan3::standard_layout_concept
 *
 * \details
 * This alphabet is not usually used directly, but instead via seqan3::masked.
 * For more information see the \link mask Mask submodule \endlink.
 *
 * \snippet test/snippet/alphabet/mask/mask.cpp general
 */
class mask : public alphabet_base<mask, 2, void>
{
private:
    //!\brief The base class.
    using base_t = alphabet_base<mask, 2, void>;

    //!\brief Befriend seqan3::alphabet_base.
    friend base_t;

public:
    /*!\name Constructors, destructor and assignment
     * \{
     */
    constexpr mask() : base_t{} {}
    constexpr mask(mask const &) = default;
    constexpr mask(mask &&) = default;
    constexpr mask & operator=(mask const &) = default;
    constexpr mask & operator=(mask &&) = default;
    ~mask() = default;
    //!\}

    /*!\name Boolean values
     * \brief Static member "booleans" that can be assigned to the alphabet or used in aggregate initialization.
     * \details Similar to an Enum interface.
     */
    //!\{
    static const mask UNMASKED;
    static const mask MASKED;
    //!\}
};

mask constexpr mask::UNMASKED{mask{}.assign_rank(0)};
mask constexpr mask::MASKED  {mask{}.assign_rank(1)};
} // namespace seqan3
