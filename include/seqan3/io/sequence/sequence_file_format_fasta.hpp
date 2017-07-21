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

/*!\file
 * \brief Contains sequence_file_in class and the tag dispatching of read function.
 * \author Temesgen H. Dadi <temesgen.dadi@fu-berlin.de>
 * \author Svenja Mehringer <svenja.mehringer@fu-berlin.de>
 * \ingroup io/sequence
 */

#pragma once

#include <vector>
#include <string>

#include <range/v3/utility/iterator.hpp>
#include <range/v3/istream_range.hpp>

#include <seqan3/io/detail/tokenization.hpp>
#include <seqan3/io/sequence/sequence_file_format.hpp>
#include <seqan3/range/container/concept.hpp>

namespace seqan3
{

struct sequence_file_format_fasta
{
public:
    static inline std::vector<std::string> file_extensions {{".fasta"},{".fa"},{".fna"}};

    //TODO make the requirements stricter
    template <typename sequence_type,
              typename meta_type,
              typename stream_type,
              typename options_type>
        requires sequence_concept<std::decay_t<sequence_type>> &&
               sequence_concept<std::decay_t<meta_type>>
    void read(sequence_type && seq,
              meta_type && meta,
              stream_type & stream,
              options_type const & options);

    template <typename seqs_type,
              typename metas_type,
              typename stream_type,
              typename options_type>
        requires sequence_of_sequence_concept<std::decay_t<seqs_type>> &&
                 sequence_of_sequence_concept<std::decay_t<metas_type>>
    void read(seqs_type && seqs,
              metas_type && metas,
              stream_type & stream,
              options_type const & options,
              size_t max_records = 0);

    template <typename sequence_type,
              typename meta_type,
              typename stream_type,
              typename options_type>
        requires sequence_concept<std::decay_t<sequence_type>> &&
                 sequence_concept<std::decay_t<meta_type>>
    void write(sequence_type && seq,
              meta_type && meta,
              stream_type & stream,
              options_type const & options);

    template <typename seqs_type,
              typename metas_type,
              typename stream_type,
              typename options_type>
              requires sequence_of_sequence_concept<std::decay_t<seqs_type>> &&
                       sequence_of_sequence_concept<std::decay_t<metas_type>>
    void write(seqs_type && seqs,
               metas_type && metas,
               stream_type & stream,
               options_type const & options,
               size_t max_records = 0);
};

template <typename sequence_type,
          typename meta_type,
          typename stream_type,
          typename options_type>
     requires sequence_concept<std::decay_t<sequence_type>> &&
              sequence_concept<std::decay_t<meta_type>>
void sequence_file_format_fasta::read(sequence_type && seq,
                                      meta_type && meta,
                                      stream_type & stream,
                                      options_type const & options)
{
    using namespace seqan3::detail;
    auto [stream_begin, stream_end] = make_preferred_input_iterator_range(stream);

    // check if we are at > and consume it
    assert(*stream_begin == '>');
    read_one(stream_begin, stream_end, std::ignore);

    // read sequence identifier
    auto meta_iter = make_preferred_output_iterator(meta);
    read_line(stream_begin, stream_end, meta_iter);

    // read sequence
    auto seq_iter = make_preferred_output_iterator(seq);
    read_until(stream_begin, stream_end, seq_iter, is_whitespace());
    read_line(stream_begin, stream_end, std::ignore);
}

template <typename seqs_type,
          typename metas_type,
          typename stream_type,
          typename options_type>
    requires sequence_of_sequence_concept<std::decay_t<seqs_type>> &&
             sequence_of_sequence_concept<std::decay_t<metas_type>>
void sequence_file_format_fasta::read(seqs_type && seqs,
                                      metas_type && metas,
                                      stream_type & stream,
                                      options_type const & options,
                                      size_t max_records)
{
    using namespace seqan3::detail;
    auto [stream_begin, stream_end] = make_preferred_input_iterator_range(stream);

    if (seqs.size() < max_records)
        seqs.resize(max_records);
    if (metas.size() < max_records)
        metas.resize(max_records);

    for (size_t i = 0; ( i < max_records && stream_begin != stream_end); ++i)
    {
        // check if we are at > and consume it
        assert(*stream_begin == '>');
        read_one(stream_begin, stream_end, std::ignore);

        // read sequence identifier
        auto meta_iter = make_preferred_output_iterator(metas[i]);
        read_line(stream_begin, stream_end, meta_iter);

        // read sequence
        auto seq_iter = make_preferred_output_iterator(seqs[i]);
        read_until(stream_begin, stream_end, seq_iter, is_whitespace());
        read_line(stream_begin, stream_end, std::ignore);
    }
}
} // namespace seqan3
