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

#include "alphabet.hpp"
#include <tuple>
#include "nucleotide/dna5.hpp"

namespace seqan3
{

template <typename first_alphabet_type, typename ...alphabet_types>
    requires alphabet_concept<first_alphabet_type> && (alphabet_concept<alphabet_types> && ...)
struct union_alphabet
{
    /* types */
    using char_type = typename first_alphabet_type::char_type;
    using integral_type = typename first_alphabet_type::integral_type;

    using tuple_type = std::tuple<first_alphabet_type, alphabet_types...>;

    template <size_t ...I>
    constexpr const integral_type* get_sizes(std::index_sequence<I...>) {
        constexpr integral_type sizes[] = {std::get<I>(tuple_type)::value_size...};
        return sizes;
    };

    integral_type value;

    constexpr char_type to_char() const
    {
        //TODO
    }

    constexpr integral_type to_integral() const
    {
        return value;
    }

    constexpr union_alphabet<first_alphabet_type, alphabet_types...> from_char(char_type const c)
    {
        //TODO
    }

    constexpr union_alphabet<first_alphabet_type, alphabet_types...> from_integral(integral_type const i)
    {
        value = i;
        return *this;
    }

    // todo autodetect type
    static constexpr integral_type value_size = first_alphabet_type::value_size + (alphabet_types::value_size + ...);
};

using alp_t = union_alphabet<dna5, dna5>;

constexpr auto alp = std::make_tuple(dna5::A, dna5::C);

static_assert(detail::internal_alphabet_concept<alp_t>);
//static_assert(alphabet_concept<alp_t>);
static_assert(all_types_are_pod(alp));

} // namespace seqan3
