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
// Author: Svenja Mehringer <svenja.mehringer@fu-berlin.de>
// Author: Temesgen H. Dadi <temesgen.dadi@fu-berlin.de>
// ==========================================================================

#pragma once

#include <vector>
#include <string>

namespace seqan3
{

struct sequence_file_format_fasta
{
public:
    static inline std::vector<std::string> file_extensions {{".fasta"},{".fa"}};

    //TODO make the requirements stricter
    template <typename sequence_type,
              typename meta_type,
              typename qual_type,
              typename stream_type,
              typename options_type>
        requires sequence_concept<std::decay_t<sequence_type>> &&
                 sequence_concept<std::decay_t<meta_type>> &&
                 sequence_concept<std::decay_t<qual_type>>
    void read(sequence_type && seq,
              meta_type && meta,
              qual_type && qual,
              stream_type & stream,
              options_type const & options);

    template <typename seqs_type,
              typename metas_type,
              typename quals_type,
              typename stream_type,
              typename options_type>
        requires sequence_of_sequence_concept<std::decay_t<seqs_type>> &&
                 sequence_of_sequence_concept<std::decay_t<metas_type>> &&
                 sequence_of_sequence_concept<std::decay_t<quals_type>>
    void read(seqs_type && seqs,
              metas_type && metas,
              quals_type && quals,
              stream_type & stream,
              options_type const & options,
              size_t max_records = 0);

    template <typename sequence_type,
              typename meta_type,
              typename qual_type,
              typename stream_type,
              typename options_type>
        requires sequence_concept<std::decay_t<sequence_type>> &&
                 sequence_concept<std::decay_t<meta_type>> &&
                 sequence_concept<std::decay_t<qual_type>>
    void write(sequence_type && seq,
               meta_type && meta,
               qual_type && qual,
               stream_type & stream,
               options_type const & options);

    template <typename seqs_type,
              typename metas_type,
              typename quals_type,
              typename stream_type,
              typename options_type>
        requires sequence_of_sequence_concept<std::decay_t<seqs_type>> &&
                 sequence_of_sequence_concept<std::decay_t<metas_type>> &&
                 sequence_of_sequence_concept<std::decay_t<quals_type>>
    void write(seqs_type && seqs,
               metas_type && metas,
               quals_type && quals,
               stream_type & stream,
               options_type const & options,
               size_t max_records = 0);
};

template <typename sequence_type,
          typename meta_type,
          typename qual_type,
          typename stream_type,
          typename options_type>
     requires sequence_concept<std::decay_t<sequence_type>> &&
              sequence_concept<std::decay_t<meta_type>> &&
              sequence_concept<std::decay_t<qual_type>>
inline void
sequence_file_format_fasta::read(sequence_type && seq,
                                 meta_type && meta,
                                 qual_type && qual,
                                 stream_type & stream,
                                 options_type const & options)
{
    // TODO:: tokenization
}

template <typename seqs_type,
          typename metas_type,
          typename quals_type,
          typename stream_type,
          typename options_type>
    requires sequence_of_sequence_concept<std::decay_t<seqs_type>> &&
             sequence_of_sequence_concept<std::decay_t<metas_type>> &&
             sequence_of_sequence_concept<std::decay_t<quals_type>>
inline void
sequence_file_format_fasta::read(seqs_type && seqs,
                                 metas_type && metas,
                                 quals_type && quals,
                                 stream_type & stream,
                                 options_type const & options,
                                 size_t max_records)
{
    // TODO:: fill containers using read()
}
} // namespace seqan3
