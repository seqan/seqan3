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
// Author: Svenja, Temesgen, Jongkyu <our_emails@fu-berlin.de>
// ==========================================================================
// sequence_file_format_fasta.hpp
// ==========================================================================
#pragma once

#include <limits>
#include <vector>
#include <string>
#include <tuple>

namespace seqan3
{

struct sequence_file_format_fasta
{
public:
    // sequence_file_format_fasta(std::string file_name)
    // {
    //     stream.open(_file_name, std::ios::in);
    // }

    // rule of six constructors deleted
    // sequence_file_format_fasta() = default;
    // sequence_file_format_fasta(sequence_file_format_fasta const & ) = default;
    // sequence_file_format_fasta & operator=(sequence_file_format_fasta const &) = default;
    // sequence_file_format_fasta(sequence_file_format_fasta &&) = default;
    // sequence_file_format_fasta & operator=(sequence_file_format_fasta &&) = default;

    static inline std::vector<std::string> file_extensions {{"fasta"},{"fa"}};

    //TODO make the requirements stricter
    template <typename sequence_type,
              typename id_type,
              typename qual_type,
              typename stream_type,
              typename options_type>
        requires container_concept<sequence_type> &&
                 container_concept<id_type> &&
                 container_concept<qual_type>
    void read(sequence_type && seq,
              id_type && id,
              qual_type && qual,
              stream_type & stream,
              options_type const & options);

    template <typename seqs_type,
              typename ids_type,
              typename quals_type,
              typename stream_type,
              typename options_type,
              size_t max_records = 0>
        requires container_concept<typename seqs_type::value> &&
                 container_concept<typename ids_type::value> &&
                 container_concept<typename quals_type::value>
    void read(seqs_type && seqs,
              ids_type && ids,
              quals_type && quals,
              stream_type & stream,
              options_type const & options);
};

/** implementations **/
template <typename sequence_type,
          typename id_type,
          typename qual_type,
          typename stream_type,
          typename options_type>
requires container_concept<sequence_type> &&
         container_concept<id_type> &&
         container_concept<qual_type>
void sequence_file_format_fasta::read(sequence_type && seq,
                                      id_type && id,
                                      qual_type && qual,
                                      stream_type & stream,
                                      options_type const & options)
{
    skip_until(stream, '>');    // forward to the next '>'
    skip_one(stream, '>');    // assert and skip '>'
    // read_line(id, options.id_filter, stream);                        // read Fasta id
    // read_until(seq, options.sequence_filter, stream, '>'); // read Fasta sequence
}
template <typename seqs_type,
          typename ids_type,
          typename quals_type,
          typename stream_type,
          typename options_type,
          size_t max_records = 0>
requires container_concept<typename seqs_type::value> &&
         container_concept<typename ids_type::value> &&
         container_concept<typename quals_type::value>
void sequence_file_format_fasta::read(seqs_type && seqs,
                                      ids_type && ids,
                                      quals_type && quals,
                                      stream_type & stream,
                                      options_type const & options)
{
    typedef typename seqs_type::value_type seq_type;
    typedef typename ids_type::value_type id_type;
    typedef typename quals_type::value_type qual_type;
    uint32_t    records_to_read = max_records;
    if (max_records == 0)
    {
        std::numeric_limits<uint32_t>::max();
    }
    for (; !stream.at_end() && records_to_read > 0; --records_to_read)
    {
        // read(seqs[last], temp_id, temp_qual, stream, options);
        // seqs.push_back(new_seq);
        // ids.push_back(new_id);
        // quals.push_back(new_qual);
    }
}
} // namespace seqan3