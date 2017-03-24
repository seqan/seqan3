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
//! The fasta Format
/*! This struct satisfies the `sequence_file_format_concept` and implements
 * read and write function for the
 * [fasta](https://www.genomatix.de/online_help/help/sequence_formats.html#FASTA) format.
 */
struct sequence_file_format_fasta
{
public:
    //! static member variable storing the file extensions
    /*!
     * The file extensions identify the format and are used for automated
     * format deduction when instantiating a `sequence_file_in` object.
     * Note: The extensions must be passed with a dot in front (e.g. ".fa")
     */
    static inline std::vector<std::string> file_extensions {{".fasta"},{".fa"}};

    //TODO make the requirements stricter
    //! reads a single record from the stream and into the given arguments
    /*!
     * \param seq the raw sequence information.
     * \param meta the meta information (e.g. the sequence identifier/name).
     * \param qual the quality information.
     * \param stream a stream to read from (e.g. std::ifstream).
     * \param options a struct containing additional options. See `seq_file_in::options_type`.
     */
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

    //! reads many or all information from the stream appending it to the given arguments
    /*!
     * \param seqs a container of sequences to append to.
     * \param metas a container of meta information to append to.
     * \param quals a container of quality information.
     * \param stream a stream to read from (e.g. std::ifstream).
     * \param options a struct containing additional options. See `seq_file_in::options_type`.
     * \param max_records limit the number of records to read to max_records.
     */
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

   //! writes a single record to the stream
   /*!
    * \param seq a container of sequences to append to.
    * \param meta a container of meta information to append to.
    * \param qual a container of quality information.
    * \param stream a stream to write to (e.g. std::ofstream).
    * \param options a struct containing additional options. See `seq_file_out::options_type`.
    */
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

   //! writes a many or all single record to the stream
   /*!
    * \param seqs a container of sequences to append to.
    * \param metas a container of meta information to append to.
    * \param quals a container of quality information.
    * \param stream a stream to write to (e.g. std::ofstream).
    * \param options a struct containing additional options. See `seq_file_out::options_type`.
    * \param max_records the number of records to write to the stream.
    */
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
