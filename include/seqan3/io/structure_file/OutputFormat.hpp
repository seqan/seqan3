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
 * \brief Provides seqan3::StructureFileOutputFormat and auxiliary classes.
 * \author Jörg Winkler <j.winkler AT fu-berlin.de>
 */

#pragma once

#include <fstream>
#include <string>
#include <vector>
#include <set>

#include <seqan3/alphabet/nucleotide/rna5.hpp>
#include <seqan3/alphabet/structure/all.hpp>
#include <seqan3/core/type_list.hpp>
#include <seqan3/io/structure_file/output_options.hpp>

namespace seqan3
{

/*!\interface seqan3::StructureFileOutputFormat <>
 * \brief The generic concept for sequence file out formats.
 * \ingroup structure_file
 *
 * \details
 *
 * The details of this concept are only relevant to developers who wish to implement their own format.
 * The requirements for this concept are given as related functions and metafunctions.
 * Types that satisfy this concept are shown as "implementing this interface".
 */
//!\cond
template <typename t>
concept StructureFileOutputFormat = requires(t & v,
                                                             std::ofstream & f,
                                                             structure_file_output_options & options,
                                                             rna5_vector & seq,
                                                             std::string & id,
                                                             std::vector<std::set<std::pair<double, size_t>>> & bpp,
                                                             std::vector<wuss51> & structure,
                                                             std::vector<structured_rna<rna5, wuss51>> & structured_seq,
                                                             double energy,
                                                             double react,
                                                             double react_err,
                                                             std::string & comment,
                                                             size_t offset)
{
    t::file_extensions;

    { v.write(f, options, seq,            id,          bpp,         structure,
                          energy,         react,       react_err,   comment,        offset)      } -> void;
    { v.write(f, options, seq,            id,          bpp,         std::ignore,
                          std::ignore,    std::ignore, std::ignore, std::ignore,    std::ignore) } -> void;
    { v.write(f, options, structured_seq, id,          std::ignore, structured_seq,
                          energy,         std::ignore, std::ignore, std::ignore,    std::ignore) } -> void;
    { v.write(f, options, std::ignore,    std::ignore, std::ignore, std::ignore,
                          std::ignore,    std::ignore, std::ignore, std::ignore,    std::ignore) } -> void;
    // the last is required to be compile time valid, but should always throw at run-time.
};
//!\endcond

/*!\name Requirements for seqan3::StructureFileOutputFormat
 * \brief You can expect these **members** on all types that implement seqan3::StructureFileOutputFormat.
 * \memberof seqan3::StructureFileOutputFormat
 * \{
 */

/*!\fn void write(stream_type & stream,
 *                structure_file_output_options const & options,
 *                seq_type && seq,
 *                id_type && id,
 *                bpp_type && bpp,
 *                structure_type && structure,
 *                energy_type && energy,
 *                react_type && react,
 *                react_type && react_err,
 *                comment_type && comment,
 *                offset_type && offset)
 * \brief Write the given fields to the specified stream.
 * \memberof seqan3::StructureFileOutputFormat
 * \tparam stream_type      Output stream, must satisfy seqan3::Ostream with `char`.
 * \tparam seq_type         Type of the seqan3::field::SEQ output; must satisfy std::ranges::OutputRange
 * over a seqan3::Alphabet.
 * \tparam id_type          Type of the seqan3::field::ID output; must satisfy std::ranges::OutputRange
 * over a seqan3::Alphabet.
 * \tparam bpp_type         Type of the seqan3::field::BPP output; must satisfy std::ranges::OutputRange
 * over a set of pair of types satisfying std::is_floating_point and std::numeric_limits::is_integer, respectively.
 * \tparam structure_type   Type of the seqan3::field::STRUCTURE output; must satisfy std::ranges::OutputRange
 * over a seqan3::RnaStructure.
 * \tparam energy_type      Type of the seqan3::field::ENERGY output; must satisfy std::is_floating_point.
 * \tparam react_type       Type of the seqan3::field::REACT and seqan3::field::REACT_ERR output;
 * must satisfy std::is_floating_point.
 * \tparam comment_type     Type of the seqan3::field::COMMENT output; must satisfy std::ranges::OutputRange
 * over a seqan3::Alphabet.
 * \tparam offset_type      Type of the seqan3::field::OFFSET output; must satisfy std::numeric_limits::is_integer.
 * \param[in,out] stream    The output stream to write into.
 * \param[in]     options   File specific options passed to the format.
 * \param[in]     seq       The data for seqan3::field::SEQ output, i.e. the "sequence".
 * \param[in]     id        The data for seqan3::field::ID output, e.g. the header line.
 * \param[in]     bpp       The data for seqan3::field::BPP output.
 * \param[in]     structure The data for seqan3::field::STRUCTURE output.
 * \param[in]     energy    The data for seqan3::field::ENERGY output.
 * \param[in]     react     The data for seqan3::field::REACT output.
 * \param[in]     react_err The data for seqan3::field::REACT_ERR output.
 * \param[in]     comment   The data for seqan3::field::COMMENT output.
 * \param[in]     offset    The data for seqan3::field::OFFSET output.
 *
 * \details
 *
 * ### Additional requirements
 *
 *   * The format must also accept std::ignore as parameter for any of the fields, however it shall throw an exception
 * if one of the fields required for writing the format is marked as such. [this shall be checked inside the function]
 *   * The format does not handle seqan3::field::STRUCTURED_SEQ, instead seqan3::structure_file_out splits it into
 * two views and passes it to the format as if they were separate.
 */
/*!\var static inline std::vector<std::string> seqan3::StructureFileOutputFormat::file_extensions
 * \brief The format type is required to provide a vector of all supported file extensions.
 */

//!\}

} // namespace seqan3

namespace seqan3::detail
{

/*!\brief Auxiliary value metafuncton that checks whether a type is a seqan3::type_list and all types meet
 * seqan3::StructureFileOutputFormat [default is false].
 * \ingroup core
 * \see seqan3::TypeListOfStructureFileOutputFormats
 */
template <typename t>
constexpr bool is_type_list_of_structure_file_output_formats_v = false;

/*!\brief Auxiliary value metafuncton that checks whether a type is a seqan3::type_list and all types meet
 * seqan3::StructureFileOutputFormat [overload].
 * \ingroup core
 * \see seqan3::TypeListOfStructureFileOutputFormats
 */
template <typename ... ts>
constexpr bool is_type_list_of_structure_file_output_formats_v<type_list<ts...>>
                = (StructureFileOutputFormat<ts> && ...);

/*!\brief Auxiliary concept that checks whether a type is a seqan3::type_list and all types meet
 * seqan3::StructureFileFormat.
 * \ingroup core
 * \see seqan3::is_type_list_of_structure_file_formats_v
 */
template <typename t>
concept TypeListOfStructureFileOutputFormats = is_type_list_of_structure_file_output_formats_v<t>;
} // namespace seqan3::detail
