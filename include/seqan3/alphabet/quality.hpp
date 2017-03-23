// ============================================================================
//                 SeqAn - The Library for Sequence Analysis
// ============================================================================
//
// Copyright (c) 2006-2017, Knut Reinert & Freie Universitaet Berlin
// Copyright (c) 2016-2017, Knut Reinert & MPI Molekulare Genetik
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
// Author: Marie Hoffmann <marie.hoffmann AT fu-berlin.de>
// ============================================================================

#pragma once

#include <iostream>
#include <string>

#include "alphabet.hpp"

namespace seqan3
{

namespace detail
{

template <typename q>
concept bool internal_quality_concept = requires (q quality)
{
    requires internal_alphabet_concept<q>;
    typename q::phred_type;

    { quality.to_phred() } -> typename q::phred_type;

    { quality.from_phred(0) } -> q;

};

} // namespace seqan3::detail

template <typename alphabet_type>
    requires detail::internal_quality_concept<alphabet_type>
struct underlying_phred
{
    using type = typename alphabet_type::phred_type;
};

template <typename alphabet_type>
using underlying_phred_t = typename underlying_phred<alphabet_type>::type;

template <typename alphabet_type>
requires detail::internal_quality_concept<alphabet_type>
constexpr alphabet_type from_phred(alphabet_type & c, char const in)
{
    return c.from_phred(in);
}

template <typename alphabet_type>
requires detail::internal_quality_concept<alphabet_type>
constexpr alphabet_char_t<alphabet_type> to_phred(alphabet_type const & c)
{
    return c.to_phred();
}

// ------------------------------------------------------------------
// concept
// ------------------------------------------------------------------

template<typename q>
concept bool quality_concept = requires(q quality)
{
    // requires fulfillment of alphabet concept
    requires alphabet_concept<q>;

    // conversion functions from t -> c and vice versa
    { from_phred(quality, typename q::integral_type{}) } -> q;
    { to_phred(quality) } -> const typename q::phred_type;
    typename underlying_phred<q>::type;
};

}
