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

#include <vector>
#include <string>

#include <seqan3/core/platform.hpp>
#include <seqan3/range/container/concept.hpp>
#include <seqan3/range/view/char_to.hpp>

namespace seqan3
{

class sequence_file_format_fasta
{
public:
    sequence_file_format_fasta() = default;
    sequence_file_format_fasta(sequence_file_format_fasta const &) = delete;
    sequence_file_format_fasta & operator=(sequence_file_format_fasta const &) = delete;
    sequence_file_format_fasta(sequence_file_format_fasta &&) = default;
    sequence_file_format_fasta & operator=(sequence_file_format_fasta &&) = default;

    static inline std::vector<std::string> file_extensions // TODO: array of constexpr_string
    {
        { "fasta" },
        { "fa" }
        // ...
    };

    template <typename record_type,
              typename stream_type, // TODO: istream_concept
              typename options_type> // TODO constrain this?
    void read(record_type & record,
              stream_type & stream,
              options_type const & /*options*/)
    {
        //TODO actual implementation goes here
        std::getline(stream, get<field::ID>(record));

        std::string buffer;
        std::getline(stream, buffer);

        get<field::SEQ>(record) = buffer | view::char_to<ranges::range_value_t<tuple_element_t<field::SEQ, record_type>>>;
    }

//     template <typename seqs_type,
//               typename ids_type,
//               typename quals_type,
//               typename stream_type,
//               typename options_type
//               size_t max_records = 0>
//         requires container_concept<typename seqs_type::value> &&
//                  container_concept<typename ids_type::value> &&
//                  container_concept<typename quals_type::value>
//     void read(seqs_type && seqs,
//               ids_type && ids,
//               quals_type && quals,
//               stream_type & stream,
//               options_type const & options);
};

/** implementations **/


} // namespace seqan3

