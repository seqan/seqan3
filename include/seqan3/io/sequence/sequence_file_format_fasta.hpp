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
/*!\brief defines FASTA format specific read and write functions
 * \ingroup sequence
 *
 * \details
 * FASTA format
 * ------------
 * FASTA format is a plain text format used mainly to sore a sequence of nucleotides or amino acids.
 * A FASTA formated file can have one or more entries. Each entry represents a contiguous sequence and is composed
 * of a sequence identifier (id) and the sequence itself. The first line of an entry starts with '>' character.
 * The rest of this line is the id of the sequence that follows. The actual sequence is written starting from a newline
 * and continues until another line starting with '>' appears or until the end of the file. Unlike the id a single
 * sequence could span multiple lines. (These lines are usually of the same length with an exception to the last line).
 *
 * ### Example #######
 *
 * > &gt;gi|186681228|ref|YP_001864424.1| phycoerythrobilin:ferredoxin oxidoreductase<br>
 * > MNSERSDVTLYQPFLDYAIAYMRSRLDLEPYPIPTGFESNSAVVGKGKNQEEVVTTSYAFQTAKLRQIRA<br>
 * > AHVQGGNSLQVLNFVIFPHLNYDLPFFGADLVTLPGGHLIALDMQPLFRDDSAYQAKYTEPILPIFHAHQ<br>
 * > QHLSWGGDFPEEAQPFFSPAFLWTRPQETAVVETQVFAAFKDYLKAYLDFVEQAEAVTDSQNLVAIKQAQ<br>
 * > LRYLRYRAEKDPARGMFKRFYGAEWTEEYIHGFLFDLERKLTVVK
 *
 *
 */
struct sequence_file_format_fasta
{
public:
    static inline std::vector<std::string> file_extensions {{".fasta"},{".fa"},{".fna"}};

    /*!\name read functions
     * \{
     */
     /*!\brief reads the next FASTA record from a file.
     * \param[out] seq a container where the fetched sequence will be written to.
     * \param[out] meta a container where the fetched sequence-id will be written to.
     * \param[in] stream an input stream to read from.
     * \param[in] options different options to be applied while reading. E.g. filters on the sequence(-id) etc.
     *
     */
    template <typename sequence_type,
              typename meta_type,
              typename stream_type,
              typename options_type>
    //!\cond
        requires sequence_concept<std::decay_t<sequence_type>> &&
               sequence_concept<std::decay_t<meta_type>>
    //!\endcond
    void read(sequence_type && seq,
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

     /*!\brief reads the next 'n' FASTA records from a file.
     * \param[out] seqs a container of containers where the fetched sequences will be written to.
     * \param[out] metas a container of containers where the fetched sequence-ids will be written to.
     * \param[in] stream an input stream to read from.
     * \param[in] options different options to be applied while reading. E.g. filters on the sequence(-id) etc.
     * \param[in] max_records the number of records to read.
     *
     */
    template <typename seqs_type,
              typename metas_type,
              typename stream_type,
              typename options_type>
    //!\cond
        requires sequence_of_sequence_concept<std::decay_t<seqs_type>> &&
                 sequence_of_sequence_concept<std::decay_t<metas_type>>
    //!\endcond
    void read(seqs_type && seqs,
              metas_type && metas,
              stream_type & stream,
              options_type const & options,
              size_t max_records = 0)
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
    //!\}

    /*!\name write functions
     * \{
     */
     /*!\brief writes a FASTA record to a file.
     * \param[in] seq a container with a sequence to be written
     * \param[in] meta a container with a sequence-id to be written
     * \param[out] stream an output stream to write to
     * \param[in] options different options to be applied before writing. E.g. line wrapping of the sequence etc.
     *
     */
    template <typename sequence_type,
              typename meta_type,
              typename stream_type,
              typename options_type>
    //!\cond
        requires sequence_concept<std::decay_t<sequence_type>> &&
                 sequence_concept<std::decay_t<meta_type>>
    //!\endcond
    void write(sequence_type && seq,
              meta_type && meta,
              stream_type & stream,
              options_type const & options);

     /*!\brief writes a batch of 'n' FASTA records to a file.
     * \param[in] seqs a container of containers with sequences to be written
     * \param[in] metas a container of containers with sequence-ids to be written
     * \param[out] stream an output stream to write to
     * \param[in] options different options to be applied before writing. E.g. line wrapping of the sequence etc.
     * \param[in] max_records the number of records to read.
     *
     */
    template <typename seqs_type,
              typename metas_type,
              typename stream_type,
              typename options_type>
    //!\cond
              requires sequence_of_sequence_concept<std::decay_t<seqs_type>> &&
                       sequence_of_sequence_concept<std::decay_t<metas_type>>
    //!\endcond
    void write(seqs_type && seqs,
               metas_type && metas,
               stream_type & stream,
               options_type const & options,
               size_t max_records = 0);
    //!\}
};

} // namespace seqan3
