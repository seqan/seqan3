// ==========================================================================
//                 SeqAn - The Library for Sequence Analysis
// ==========================================================================
//
// Copyright (c) 2006-2017, Knut Reinert, FU Berlin
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
// ==========================================================================
// Author: Hannes Hauswedell <hannes.hauswedell@fu-berlin.de>
// Author: Marcel Ehrhardt <marcel.ehrhardt@fu-berlin.de>
// Author: David Heller <david.heller@fu-berlin.de>
// ==========================================================================

#pragma once

#include "gap.hpp"
#include "../union_alphabet.hpp"
#include "../nucleotide/dna4.hpp"

namespace seqan3
{

template <typename underlying_t>
    requires alphabet_concept<underlying_t>
struct gapped_alphabet : public union_alphabet<underlying_t, gap>
{
    using union_alphabet<underlying_t, gap>::_value;
    using union_alphabet<underlying_t, gap>::value_size;

    using union_alphabet<underlying_t, gap>::union_alphabet;
    using union_alphabet<underlying_t, gap>::operator=;

    using typename union_alphabet<underlying_t, gap>::rank_type;
    using typename union_alphabet<underlying_t, gap>::char_type;

    /* public member functions */
    constexpr bool is_gap() const
    {
        return _value == value_size - 1;
    }

    constexpr gapped_alphabet set_gap()
    {
        _value = value_size - 1;
        return *this;
    }

    constexpr gapped_alphabet & assign_rank(rank_type const i)
    {
        union_alphabet<underlying_t, gap>::assign_rank(i);
        return *this;
    }

    constexpr gapped_alphabet & assign_char(char_type const c)
    {
        union_alphabet<underlying_t, gap>::assign_char(c);
        return *this;
    }
};

#ifndef NDEBUG
static_assert(alphabet_concept<gapped_alphabet<dna4>>);
#endif

}
