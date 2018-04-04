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
 * \brief Provides the seqan3::sequence_file_format_fasta class.
 * \author Hannes Hauswedell <hannes.hauswedell AT fu-berlin.de>
 */

#pragma once

#include <string>
#include <vector>

#include <seqan3/alphabet/nucleotide/dna5.hpp>
#include <seqan3/alphabet/quality/aliases.hpp>
#include <seqan3/range/concept.hpp>
#include <seqan3/range/view/char_to.hpp>

namespace seqan3
{
// will be documented when properly implemented
//!\cond
class sequence_file_format_fasta
{
public:
    sequence_file_format_fasta() = default;
    sequence_file_format_fasta(sequence_file_format_fasta const &) = delete;
    sequence_file_format_fasta & operator=(sequence_file_format_fasta const &) = delete;
    sequence_file_format_fasta(sequence_file_format_fasta &&) = default;
    sequence_file_format_fasta & operator=(sequence_file_format_fasta &&) = default;

    static inline std::vector<std::string> file_extensions
    {
        { "fasta" },
        { "fa" }
        // ...
    };

    template <typename stream_type,     // constraints checked by file
              typename options_type,    // given by file
              typename seq_type,        // other constraints checked inside function
              typename id_type,
              typename qual_type,
              typename seq_qual_type>
    void read(stream_type        & stream,
              [[maybe_unused]] options_type const & options,
              [[maybe_unused]] seq_type           & seq,
              [[maybe_unused]] id_type            & id,
              [[maybe_unused]] qual_type          & qual,
              [[maybe_unused]] seq_qual_type      & seq_qual)
    {
        //TODO check if seq_qual or seq+qual
        //TODO create output iterators on ranges
        //TODO pass to actual implementation

        std::string buffer;
        std::getline(stream, buffer);

        if constexpr (output_range_concept<id_type, char>)
            id = buffer;

        std::getline(stream, buffer);

        if constexpr (output_range_concept<seq_type, dna5>)
            seq = buffer | view::char_to<dna5>;

        if constexpr (output_range_concept<seq_qual_type, dna5q>)
            seq_qual = buffer | view::char_to<dna5q>;
    }

    template <typename stream_type,     // constraints checked by file
              typename options_type,    // given by file
              typename seq_type,        // other constraints checked inside function
              typename id_type,
              typename qual_type,
              typename seq_qual_type>
    void write([[maybe_unused]] stream_type        & stream,
               [[maybe_unused]] options_type const & options,
               [[maybe_unused]] seq_type           && seq,
               [[maybe_unused]] id_type            && id,
               [[maybe_unused]] qual_type          && qual,
               [[maybe_unused]] seq_qual_type      && seq_qual)
    {
        //TODO
    }
};
//!\endcond

} // namespace seqan3

