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

// ==================================================================
// sequence_file_format_fasta
// ==================================================================

//! The fasta file format.
/*! This struct satisfies the `sequence_file_format_concept` and implements
 * read and write function for the
 * [fasta](https://www.genomatix.de/online_help/help/sequence_formats.html#FASTA) format.
 */
struct sequence_file_format_fasta
{
public:
    //! Static member variable storing the valid file extensions.
    /*!
     * The valid file extensions that are used to identify the fasta format.
     * Note: The extensions must be passed with a dot in front (e.g. ".fa")
     */
    static inline std::vector<std::string> file_extensions {{".fasta"}, {".fa"}};

    //TODO make the requirements stricter
    //! Reads a single fasta record from the stream.
    /*!
     * \param[in,out] seq The raw sequence information.
     * \param[in,out] meta The meta information (e.g. the sequence identifier/name).
     * \param[in,out] qual The quality information.
     * \param[in,out] stream The input stream to read from (e.g. std::ifstream).
     * \param[in]     options A struct containing additional options. See `seq_file_in::options_type`.
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

    //! Reads a chunk or all fasta records from the stream and appends it to the given container.
    /*!
     * \param[in,out] seqs A container of sequences to append to.
     * \param[in,out] metas A container of meta information to append to.
     * \param[in,out] quals A container of quality information.
     * \param[in,out] stream The input stream to read from (e.g. std::ifstream).
     * \param[in]     options A struct containing additional options. See `seq_file_in::options_type`.
     * \param[in]     max_records Limits the number of records to read to max_records. Defaults to `0` 
     *                            indicating that all records should be read.
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
              size_t const max_records = 0);

   //! Writes a single fasta record to the stream.
   /*!
    * \param[in,out] The output stream to write to (e.g. std::ofstream).
    * \param[in] seq A container of sequences to append to.
    * \param[in] meta A container of meta information to append to.
    * \param[in] qual A container of quality information.
    * \param[in] options A struct containing additional options. See `seq_file_out::options_type`.
    */
    template <typename stream_type,
              typename sequence_type,
              typename meta_type,
              typename qual_type,
              typename options_type>
        requires sequence_concept<std::decay_t<sequence_type>> &&
                 sequence_concept<std::decay_t<meta_type>> &&
                 sequence_concept<std::decay_t<qual_type>>
    void write(stream_type & stream,
               sequence_type const & seq,
               meta_type const & meta,
               qual_type const & qual,
               options_type const & options);

   //! Writes a chunk or all fasta records to the stream.
   /*!
    * \param[in,out] The output stream to write to (e.g. std::ofstream).
    * \param[in] seqs A container of sequences to append to.
    * \param[in] metas A container of meta information to append to.
    * \param[in] quals A container of quality information.
    * \param[in] options A struct containing additional options. See `seq_file_out::options_type`.
    * \param[in] max_records The number of records to write to the stream. The value `0` 
    *                        indicates that all records are written. Defaults to `0`.
    */
    template <typename stream_type,
              typename seqs_type,
              typename metas_type,
              typename quals_type,
              typename options_type>
        requires sequence_of_sequence_concept<std::decay_t<seqs_type>> &&
                 sequence_of_sequence_concept<std::decay_t<metas_type>> &&
                 sequence_of_sequence_concept<std::decay_t<quals_type>>
    void write(stream_type & stream,
               seqs_type const & seqs,
               metas_type const & metas,
               quals_type const & quals,
               options_type const & options,
               size_t const max_records = 0);
};

// ==================================================================
// public functions
// ==================================================================

// ------------------------------------------------------------------
// function read
// ------------------------------------------------------------------

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
