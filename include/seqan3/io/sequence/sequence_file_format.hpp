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
// Author: Joerg Winkler <j.winkler AT fu-berlin.de>
// ============================================================================

#pragma once

#include <seqan3/alphabet/nucleotide/dna4.hpp>
#include <seqan3/io/sequence/sequence_file_in.hpp>
#include <string>
#include <fstream>

namespace seqan3
{

template <typename t>
concept bool sequence_file_format_concept = requires (t v)
{
    t::file_extensions;

    {
        v.read(dna4_string{},     // sequence
               std::string{},     // meta
               std::string{},     // quality
               std::ifstream{},   // stream
               options_type{})    // options
    };

    {
        v.write(dna4_string{},    // sequence
                std::string{},    // meta
                std::string{},    // quality
                std::ofstream{},  // stream
                options_type{})   // options
    };

};

namespace detail
{

template <size_t index, typename variant_type>
constexpr bool meets_concept_sequence_file_format()
{
    if constexpr (index == variant_size_v<variant_type>)
        return true;
    else if constexpr (!sequence_file_format_concept<variant_alternative_t<index, variant_type>>)
        return false;
    else
        return meets_concept_sequence_file_format<index+1, variant_type>();
}

} // namespace detail

} // namespace seqan3
