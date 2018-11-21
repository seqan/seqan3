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
 * \brief Provides seqan3::SequenceFileFormatOut and auxiliary classes.
 * \author Hannes Hauswedell <hannes.hauswedell AT fu-berlin.de>
 */

#pragma once

#include <fstream>
#include <string>
#include <vector>

#include <seqan3/alphabet/nucleotide/dna5.hpp>
#include <seqan3/alphabet/quality/aliases.hpp>
#include <seqan3/alphabet/quality/phred42.hpp>
#include <seqan3/core/type_list.hpp>
#include <seqan3/io/sequence_file/output_options.hpp>

namespace seqan3
{

/*!\interface seqan3::SequenceFileOutputFormat <>
 * \brief The generic concept for sequence file out formats.
 * \ingroup sequence
 *
 * \details
 *
 * The details of this concept are only relevant to developers who wish to implement their own format.
 * The requirements for this concept are given as related functions and metafunctions.
 * Types that satisfy this concept are shown as "implementing this interface".
 */
//!\cond
template <typename t>
concept SequenceFileOutputFormat = requires (t                         & v,
                                                          std::ofstream             & f,
                                                          sequence_file_output_options & options,
                                                          dna5_vector               & seq,
                                                          std::string               & id,
                                                          std::vector<phred42>   & qual,
                                                          std::vector<dna5q>        & seq_qual)
{
    t::file_extensions;

    { v.write(f, options, seq,         id,          qual)        } -> void;
    { v.write(f, options, std::ignore, id,          std::ignore) } -> void;
    { v.write(f, options, std::ignore, std::ignore, std::ignore) } -> void;
    // the last is required to be compile time valid, but should always throw at run-time.
};
//!\endcond

/*!\name Requirements for seqan3::SequenceFileOutputFormat
 * \brief You can expect these **members** on all types that implement seqan3::SequenceFileOutputFormat.
 * \memberof seqan3::SequenceFileOutputFormat
 * \{
 */

/*!\fn void write(stream_type & stream, seqan3::sequence_file_output_options const & options, seq_type && sequence,
 *                id_type && id, qual_type && qualities)
 * \brief Write the given fields to the specified stream.
 * \memberof seqan3::SequenceFileOutputFormat
 * \tparam stream_type      Output stream, must satisfy seqan3::Ostream with `char`.
 * \tparam seq_type         Type of the seqan3::field::SEQ output; must satisfy std::ranges::OutputRange
 * over a seqan3::Alphabet.
 * \tparam id_type          Type of the seqan3::field::ID output; must satisfy std::ranges::OutputRange
 * over a seqan3::Alphabet.
 * \tparam qual_type        Type of the seqan3::field::QUAL output; must satisfy std::ranges::OutputRange
 * over a seqan3::Quality.
 * \param[in,out] stream    The output stream to write into.
 * \param[in]     options   File specific options passed to the format.
 * \param[in]     sequence  The data for seqan3::field::SEQ, i.e. the "sequence".
 * \param[in]     id        The data for seqan3::field::ID, e.g. the header line in FastA.
 * \param[in]     qualities The data for seqan3::field::QUAL.
 *
 * \details
 *
 * ### Additional requirements
 *
 *   * The format must also accept std::ignore as parameter for any of the fields, however it shall throw an exception
 * if one of the fields required for writing the format is marked as such. [this shall be checked inside the function]
 *   * The format does not handle seqan3::field::SEQ_QUAL, instead seqan3::sequence_file_output splits it into two views
 *     and passes it to the format as if they were separate.
 */
/*!\var static inline std::vector<std::string> seqan3::SequenceFileOutputFormat::file_extensions
 * \brief The format type is required to provide a vector of all supported file extensions.
 */

//!\}

} // namespace seqan3

namespace seqan3::detail
{

/*!\brief Auxiliary value metafuncton that checks whether a type is a seqan3::type_list and all types meet
 * seqan3::SequenceFileOutputFormat [default is false].
 * \ingroup core
 * \see seqan3::TypeListOfSequenceFileOutputFormats
 */
template <typename t>
constexpr bool is_type_list_of_sequence_file_output_formats_v = false;

/*!\brief Auxiliary value metafuncton that checks whether a type is a seqan3::type_list and all types meet
 * seqan3::SequenceFileOutputFormat [overload].
 * \ingroup core
  * \see seqan3::TypeListOfSequenceFileOutputFormats
 */
template <typename ... ts>
constexpr bool is_type_list_of_sequence_file_output_formats_v<type_list<ts...>> =
    (SequenceFileOutputFormat<ts> && ...);

/*!\brief Auxiliary concept that checks whether a type is a seqan3::type_list and all types meet seqan3::SequenceFileFormat.
 * \ingroup core
 * \see seqan3::is_type_list_of_sequence_file_formats_v
 */
template <typename t>
concept TypeListOfSequenceFileOutputFormats = is_type_list_of_sequence_file_output_formats_v<t>;
} // namespace seqan3::detail
