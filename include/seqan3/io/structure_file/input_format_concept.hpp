// -----------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2019, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2019, Knut Reinert & MPI für molekulare Genetik
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at: https://github.com/seqan/seqan3/blob/master/LICENSE.md
// -----------------------------------------------------------------------------------------------------

/*!\file
 * \brief Provides seqan3::StructureFileInputFormat.
 * \author Jörg Winkler <j.winkler AT fu-berlin.de>
 */

#pragma once

#include <fstream>
#include <set>
#include <string>
#include <utility>
#include <vector>

#include <seqan3/alphabet/nucleotide/rna5.hpp>
#include <seqan3/alphabet/structure/structured_rna.hpp>
#include <seqan3/alphabet/structure/wuss.hpp>
#include <seqan3/core/type_list.hpp>
#include <seqan3/io/structure_file/input_options.hpp>

namespace seqan3
{
/*!\interface seqan3::StructureFileInputFormat <>
 * \brief The generic concept for structure file in formats.
 * \ingroup structure_file
 *
 * \details
 *
 * The details of this concept are only relevant to developers who wish to implement their own format.
 * The requirements for this concept are given as related functions and metafunctions.
 * Types that satisfy this concept are shown as "implementing this interface".
 */
//!\cond
template<typename t>
SEQAN3_CONCEPT StructureFileInputFormat = requires(t & v,
                                                   std::ifstream & f,
                                                   structure_file_input_options<rna5, false> & options,
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

/*!\name Requirements for seqan3::StructureFileInputFormat
 * \brief You can expect these **members** on all types that implement seqan3::StructureFileInputFormat.
 * \memberof seqan3::StructureFileInputFormat
 * \{
 */
/*!\fn void read(stream_type & stream,
 *               structure_file_input_options<seq_legal_alph_type, structured_seq_combined> const & options,
 *               seq_type & seq,
 *               id_type & id,
 *               bpp_type & bpp,
 *               structure_type & structure,
 *               energy_type & energy,
 *               react_type & react,
 *               react_type & react_err,
 *               comment_type & comment,
 *               offset_type & offset)
 * \brief Read from the specified stream and back-insert into the given field buffers.
 * \tparam stream_type      Input stream, must satisfy seqan3::istream_concept with `char`.
 * \tparam seq_type         Type of the seqan3::field::SEQ input; must satisfy std::ranges::OutputRange
 * over a seqan3::Alphabet.
 * \tparam id_type          Type of the seqan3::field::ID input; must satisfy std::ranges::OutputRange
 * over a seqan3::Alphabet.
 * \tparam bpp_type         Type of the seqan3::field::BPP input; must satisfy std::ranges::OutputRange
 * over a set of pair of types satisfying std::is_floating_point and std::numeric_limits::is_integer, respectively.
 * \tparam structure_type   Type of the seqan3::field::STRUCTURE input; must satisfy std::ranges::OutputRange
 * over a seqan3::RnaStructureAlphabet.
 * \tparam energy_type      Type of the seqan3::field::ENERGY input; must satisfy std::is_floating_point.
 * \tparam react_type       Type of the seqan3::field::REACT and seqan3::field::REACT_ERR input;
 * must satisfy std::is_floating_point.
 * \tparam comment_type     Type of the seqan3::field::COMMENT input; must satisfy std::ranges::OutputRange
 * over a seqan3::Alphabet.
 * \tparam offset_type      Type of the seqan3::field::OFFSET input; must satisfy std::numeric_limits::is_integer.
 * \param[in,out] stream    The input stream to read from.
 * \param[in]     options   File specific options passed to the format.
 * \param[out]    seq       The buffer for seqan3::field::SEQ input, i.e. the "sequence".
 * \param[out]    id        The buffer for seqan3::field::ID input, e.g. the header line.
 * \param[out]    bpp       The buffer for seqan3::field::BPP input.
 * \param[out]    structure The buffer for seqan3::field::STRUCTURE input.
 * \param[out]    energy    The buffer for seqan3::field::ENERGY input.
 * \param[out]    react     The buffer for seqan3::field::REACT input.
 * \param[out]    react_err The buffer for seqan3::field::REACT_ERR input.
 * \param[out]    comment   The buffer for seqan3::field::COMMENT input.
 * \param[out]    offset    The buffer for seqan3::field::OFFSET input.
 *
 * \details
 *
 * ### Additional requirements
 *
 *   * The function must also accept std::ignore as parameter for any of the fields.
 *     [this is enforced by the concept checker!]
 *   * In this case the data read for that field shall be discarded by the format.
 *   * Instead of passing the fields seqan3::field::SEQ and seqan3::field::STRUCTURE, you may also pass
 *     seqan3::field::STRUCTURED_SEQ to both parameters. If you do, the seqan3::value_type_t of the argument must be
 *     a specialisation of seqan3::structured_rna and the second template parameter to
 *     seqan3::structure_file_input_options must be set to true.
 */
 /*!\var static inline std::vector<std::string> seqan3::StructureFileInputFormat::file_extensions
 * \brief The format type is required to provide a vector of all supported file extensions.
 */
//!\}

} // namespace seqan3

namespace seqan3::detail
{

/*!\brief Auxiliary value metafuncton that checks whether a type is a seqan3::type_list and all types meet
 * seqan3::structure_file_format_concept [default is false].
 * \ingroup core
 * \see seqan3::TypeListOfStructureFileInputFormats
 */
template<typename t>
constexpr bool is_type_list_of_structure_file_input_formats_v = false;

/*!\brief Auxiliary value metafuncton that checks whether a type is a seqan3::type_list and all types meet
 * seqan3::StructureFileInputFormat [overload].
 * \ingroup core
 * \see seqan3::TypeListOfStructureFileInputFormats
 */
template<typename ... ts>
constexpr bool is_type_list_of_structure_file_input_formats_v<type_list<ts...>>
                = (StructureFileInputFormat<ts> && ...);

/*!\brief Auxiliary concept that checks whether a type is a seqan3::type_list and all types meet
 * seqan3::StructureFileInputFormat.
 * \ingroup core
 * \see seqan3::is_type_list_of_structure_file_formats_v
 */
template<typename t>
SEQAN3_CONCEPT TypeListOfStructureFileInputFormats = is_type_list_of_structure_file_input_formats_v<t>;
} // namespace seqan3::detail
