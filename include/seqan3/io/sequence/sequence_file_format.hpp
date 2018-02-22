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
 * \ingroup io
 * \author Joerg Winkler <j.winkler AT fu-berlin.de>
 * \brief Defines the concept of sequence file formats.
 */

#pragma once

#include <string>
#include <fstream>

#include <seqan3/alphabet/nucleotide/dna4.hpp>
#include <seqan3/io/sequence/sequence_file_in.hpp>

#if 0
//TODO(rrahn): this is a prototype and needs more refinement, disabling for now
//!\cond
namespace seqan3
{

/*!\name Requirements for seqan3::sequence_file_format_concept
 * \brief You can expect these functions on all types that implement seqan3::sequence_file_format_concept.
 * \relates seqan3::sequence_file_format_concept
 * \{
 */

/*!\fn void seqan3::sequence_file_format_concept::read(dna4_string sequence, std::string meta, std::string quality, std::ifstream stream, options_type opt)
 * \brief The read function reads a sequence file from ifstream and writes the contents into three buffers.
 * \param sequence This is the buffer for sequences. It should support any alphabet, requires at least dna4_string.
 * \param meta This is the buffer for meta-information (e.g. fasta header).
 * \param quality This is the buffer for quality values (e.g. from fastq files).
 * \param stream This is the input file stream, i.e. the source of data.
 * \param opt For passing some further options or data types to the function.
 */
/*!\fn void seqan3::sequence_file_format_concept::write(dna4_string sequence, std::string meta, std::string quality, std::ofstream stream, options_type opt)
 * \brief The write function writes the contents of sequence, meta, and quality to a file via ofstream.
 * \param sequence This is the sequence source. It should support any alphabet, requires at least dna4_string.
 * \param meta This is the source for meta-information (e.g. fasta header).
 * \param quality This is the source for quality values (e.g. from fastq files).
 * \param stream This is the output file stream, i.e. the target file.
 * \param opt For passing some further options or data types to the function.
 */
/*!\var static std::vector<std::string> seqan3::sequence_file_format_concept::file_extensions
 * \brief The format type is required to have a vector of all supported file extensions.
 */

//!\}

/*!\interface seqan3::sequence_file_format_concept <>
 * \brief The generic concept for sequence file formats.
 * \tparam t The sequence file format type.
 * \ingroup io
 *
 * The requirements for this concept are given as related functions and metafunctions.
 * Types that satisfy this concept are shown as "implementing this interface".
 */
//!\cond
template <typename t>
concept bool sequence_file_format_concept = requires (t v)
{
    { t::file_extensions; } -> std::vector<std:string>;

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
//!\endcond

namespace detail
{
//!\privatesection
/*!
 * \brief Check whether all given formats meet the sequence_file_format_concept.
 * \tparam indices Indices to access the sequence formats in the std::variant type.
 * \tparam variant_type A std::variant type that contains the format types to be checked
 * (e.g. [seqan3::sequence_file_format_fasta]).
 * \sa seqan3::detail::all_satisfy_sequence_file_format_concept()
 * \attention This is only a helper for seqan3::detail::all_satisfy_sequence_file_format_concept(),
 * please do not call it directly.
 * \returns False if at least one format does not meet the sequence file format concept, true otherwise.
 */
template<typename variant_type, std::size_t... indices>
constexpr bool all_satisfy_sequence_file_format_concept_impl(std::index_sequence<indices...>)
{
    return (sequence_file_format_concept<std::variant_alternative_t<indices, variant_type>> && ... && true);
}

/*!
 * \brief Check whether all given formats meet the sequence_file_format_concept.
 * \tparam variant_type A std::variant type that contains the format types to be checked
 * (e.g. [seqan3::sequence_file_format_fasta]).
 * \returns False if at least one format does not meet the sequence file format concept, true otherwise.
 * \ingroup io
 *
 * Usage example:
 * ```cpp
 * bool valid = detail::all_satisfy_sequence_file_format_concept<typename t::valid_formats>();
 * ```
 */
template<typename variant_type>
constexpr bool all_satisfy_sequence_file_format_concept()
{
    using index_sequence_t = std::make_integer_sequence<std::variant_size_v<variant_type>>;
    return all_satisfy_sequence_file_format_concept_impl<variant_type>(index_sequence_t{});
}
//!\publicsection

} // namespace detail

} // namespace seqan3

//!\endcond
#endif
