// -----------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2019, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2019, Knut Reinert & MPI für molekulare Genetik
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at: https://github.com/seqan/seqan3/blob/master/LICENSE.md
// -----------------------------------------------------------------------------------------------------

/*!\file
 * \brief Provides seqan3::SequenceFileInputFormat and auxiliary classes.
 * \author Jörg Winkler <j.winkler AT fu-berlin.de>
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
#include <seqan3/io/sequence_file/input_options.hpp>

namespace seqan3::detail
{

//!\brief The sequence file input format base class.
template <typename format_tag>
class sequence_file_input_format
{};

} // namespace seqan3::detail

namespace seqan3
{

/*!\interface seqan3::SequenceFileInputFormat <>
 * \brief The generic concept for sequence file in formats.
 * \ingroup sequence
 *
 * \details
 *
 * The details of this concept are only relevant to developers who wish to implement their own format.
 * The requirements for this concept are given as related functions and type traits.
 * Types that satisfy this concept are shown as "implementing this interface".
 */
//!\cond
template <typename t>
SEQAN3_CONCEPT SequenceFileInputFormat = requires (detail::sequence_file_input_format<t>    & v,
                                                   std::ifstream                            & f,
                                                   sequence_file_input_options<dna5, false> & options,
                                                   dna5_vector                              & seq,
                                                   std::string                              & id,
                                                   std::vector<phred42>                     & qual,
                                                   std::vector<dna5q>                       & seq_qual)
{
    t::file_extensions;

    { v.read(f, options, seq,         id,          qual)        } -> void;
    { v.read(f, options, seq_qual,    id,          seq_qual)    } -> void;
    { v.read(f, options, std::ignore, std::ignore, std::ignore) } -> void;
};
//!\endcond

/*!\name Requirements for seqan3::SequenceFileInputFormat
 * \brief You can expect these **members** on all types that implement seqan3::SequenceFileInputFormat.
 * \memberof seqan3::SequenceFileInputFormat
 * \{
 */

/*!\fn void read(stream_type & stream, seqan3::sequence_file_input_options const & options, seq_type & sequence,
 *               id_type & id, qual_type & qualities)
 * \brief Read from the specified stream and back-insert into the given field buffers.
 * \tparam stream_type      Input stream, must satisfy seqan3::IStream with `char`.
 * \tparam seq_type         Type of the seqan3::field::SEQ input; must satisfy std::ranges::OutputRange
 * over a seqan3::Alphabet.
 * \tparam id_type          Type of the seqan3::field::ID input; must satisfy std::ranges::OutputRange
 * over a seqan3::Alphabet.
 * \tparam qual_type        Type of the seqan3::field::QUAL input; must satisfy std::ranges::OutputRange
 * over a seqan3::WritableQualityAlphabet.
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
 *   * The function must also accept std::ignore as parameter for any of the fields.
 *     [This is enforced by the concept checker!]
 *   * In this case the data read for that field shall be discarded by the format.
 *   * Instead of passing the fields seqan3::field::SEQ and seqan3::field::QUAL, you may also pass
 *     seqan3::field::SEQ_QUAL to both parameters. If you do, the seqan3::value_type_t of the argument must be
 *     a specialisation of seqan3::qualified and the second template parameter to
 *     seqan3::sequence_file_input_options must be set to true.
 */
 /*!\var static inline std::vector<std::string> seqan3::SequenceFileInputFormat::file_extensions
 * \brief The format type is required to provide a vector of all supported file extensions.
 */
//!\}

} // namespace seqan3

namespace seqan3::detail
{

/*!\brief Auxiliary value metafuncton that checks whether a type is a seqan3::type_list and all types meet
 * seqan3::SequenceFileInputFormat [default is false].
 * \ingroup core
 * \see seqan3::TypeListOfSequenceFileInputFormats
 */
template <typename t>
constexpr bool is_type_list_of_sequence_file_input_formats_v = false;

/*!\brief Auxiliary value metafuncton that checks whether a type is a seqan3::type_list and all types meet
 * seqan3::SequenceFileInputFormat [overload].
 * \ingroup core
  * \see seqan3::TypeListOfSequenceFileInputFormats
 */
template <typename ... ts>
constexpr bool is_type_list_of_sequence_file_input_formats_v<type_list<ts...>> =
    (SequenceFileInputFormat<ts> && ...);

/*!\brief Auxiliary concept that checks whether a type is a seqan3::type_list and all types meet
 *        seqan3::SequenceFileInputFormat.
 * \ingroup core
 * \see seqan3::is_type_list_of_sequence_file_formats_v
 */
template <typename t>
SEQAN3_CONCEPT TypeListOfSequenceFileInputFormats = is_type_list_of_sequence_file_input_formats_v<t>;
} // namespace seqan3::detail
