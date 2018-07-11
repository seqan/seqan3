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
 * \brief Provides seqan3::structure_file_(in|out)_format_concept and options.
 * \author Jörg Winkler <j.winkler AT fu-berlin.de>
 */

#pragma once

#include <fstream>
#include <set>
#include <string>
#include <vector>

#include <seqan3/alphabet/nucleotide/rna5.hpp>
#include <seqan3/alphabet/structure/wuss.hpp>
#include <seqan3/core/type_list.hpp>

namespace seqan3
{
/*!\brief The options type defines various option members that influence the behaviour of all or some formats.
 * \tparam seq_legal_alphabet_ The sequence legal alphabet exposed as type trait to the format.
 * \tparam structure_legal_alphabet_ The structure legal alphabet exposed as type trait to the format.
 * \tparam structured_seq_combined Trait that exposes to the format whether seq and qual arguments are actually the
 * same/combined.
 */
template<typename seq_legal_alphabet_, bool structured_seq_combined>
struct structure_file_in_options
{
    //!\brief Export the (possibly larger) legal sequence alphabet to the format.
    using seq_legal_alphabet = seq_legal_alphabet_;

    //!\brief Read the ID string only up until the first whitespace character.
    bool truncate_ids = false;
};

/*!\interface seqan3::structure_file_in_format_concept <>
 * \brief The generic concept for structure file formats.
 * \ingroup structure
 *
 * \details
 *
 * The details of this concept are only relevant to developers who wish to implement their own format.
 * The requirements for this concept are given as related functions and metafunctions.
 * Types that satisfy this concept are shown as "implementing this interface".
 */
//!\cond
template<typename t>
concept bool structure_file_in_format_concept = requires(t & v,
                                                         std::ifstream & f,
                                                         structure_file_in_options<rna5, false> & options,
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

    { v.read(f, options, seq,            id,          bpp,         structure,
                         energy,         react,       react_err,   comment,        offset)      } -> void;

    { v.read(f, options, seq,            id,          bpp,         std::ignore,
                         std::ignore,    std::ignore, std::ignore, std::ignore,    std::ignore) } -> void;

    { v.read(f, options, structured_seq, id,          std::ignore, structured_seq,
                         energy,         std::ignore, std::ignore, std::ignore,    std::ignore) } -> void;

    { v.read(f, options, std::ignore,    std::ignore, std::ignore, std::ignore,
                         std::ignore,    std::ignore, std::ignore, std::ignore,    std::ignore) } -> void;
    // the last is required to be compile time valid, but should always throw at run-time.
};
//!\endcond

/*!\name Requirements for seqan3::structure_file_in_format_concept
 * \brief You can expect these **members** on all types that implement seqan3::structure_file_in_format_concept.
 * \memberof seqan3::structure_file_in_format_concept
 * \{
 */

// TODO(joergi-w) documentation
/*!\fn void read(stream_type & stream, seqan3::structure_file_in_options const & options, seq_type & sequence,
 *               id_type & id, qual_type & qualities)
 * \brief Read from the specified stream and back-insert into the given field buffers.
 * \memberof seqan3::sequence_file_in_format_concept
 * \tparam stream_type      Input stream, must satisfy seqan3::istream_concept with `char`.
 * \tparam seq_type         Type of the seqan3::field::SEQ input; must satisfy seqan3::output_range_concept
 * over a seqan3::alphabet_concept.
 * \tparam id_type          Type of the seqan3::field::ID input; must satisfy seqan3::output_range_concept
 * over a seqan3::alphabet_concept.
 * \tparam qual_type        Type of the seqan3::field::QUAL input; must satisfy seqan3::output_range_concept
 * over a seqan3::quality_concept.
 * \param[in,out] stream    The input stream to read from.
 * \param[in]     options   File specific options passed to the format.
 * \param[out]    sequence  The buffer for seqan3::field::SEQ input, i.e. the "sequence".
 * \param[out]    id        The buffer for seqan3::field::ID input, e.g. the header line in FastA.
 * \param[out]    qualities The buffer for seqan3::field::QUAL input.
 *
 * \details
 *
 * ### Additional requirements
 *
 *   * The function must also accept std::ignore as parameter for any of the fields. [this is enforced by the concept checker!]
 *   * In this case the data read for that field shall be discarded by the format.
 *   * Instead of passing the fields seqan3::field::SEQ and seqan3::field::QUAL, you may also pass
 *     seqan3::field::SEQ_QUAL to both parameters. If you do, the seqan3::value_type_t of the argument must be
 *     a specialisation of seqan3::quality_composition and the second template parameter to
 *     seqan3::sequence_file_in_options must be set to true.
 */
//!\}

} // namespace seqan3

namespace seqan3::detail
{

/*!\brief Auxiliary value metafuncton that checks whether a type is a seqan3::type_list and all types meet
 * seqan3::structure_file_format_concept [default is false].
 * \ingroup core
 * \see seqan3::type_list_of_structure_file_in_formats_concept
 */
template<typename t>
constexpr bool is_type_list_of_structure_file_in_formats_v = false;

/*!\brief Auxiliary value metafuncton that checks whether a type is a seqan3::type_list and all types meet
 * seqan3::structure_file_in_format_concept [overload].
 * \ingroup core
  * \see seqan3::type_list_of_structure_file_in_formats_concept
 */
template<typename ... ts>
constexpr bool is_type_list_of_structure_file_in_formats_v<type_list<ts...>>
                = (structure_file_in_format_concept<ts> && ...);

/*!\brief Auxiliary concept that checks whether a type is a seqan3::type_list and all types meet
 * seqan3::structure_file_in_format_concept.
 * \ingroup core
 * \see seqan3::is_type_list_of_structure_file_formats_v
 */
template<typename t>
concept bool type_list_of_structure_file_in_formats_concept = is_type_list_of_structure_file_in_formats_v<t>;

} // namespace seqan3::detail
