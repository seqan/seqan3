// SPDX-FileCopyrightText: 2006-2025 Knut Reinert & Freie Universität Berlin
// SPDX-FileCopyrightText: 2016-2025 Knut Reinert & MPI für molekulare Genetik
// SPDX-License-Identifier: BSD-3-Clause

/*!\file
 * \brief Provides seqan3::sequence_file_input_format and auxiliary classes.
 * \author Jörg Winkler <j.winkler AT fu-berlin.de>
 * \author Hannes Hauswedell <hannes.hauswedell AT fu-berlin.de>
 */

#pragma once

#include <fstream>
#include <string>
#include <vector>

#include <seqan3/alphabet/nucleotide/dna5.hpp>
#include <seqan3/alphabet/quality/phred42.hpp>
#include <seqan3/alphabet/quality/qualified.hpp>
#include <seqan3/io/sequence_file/input_options.hpp>
#include <seqan3/utility/type_list/type_list.hpp>

namespace seqan3::detail
{

/*!\brief Internal class used to expose the actual format interface to read sequence records from the file.
 * \ingroup io_sequence_file
 *
 * \tparam format_type The type of the format to be exposed.
 *
 * \details
 *
 * Exposes the protected member function `read_sequence_record` from the given `format_type`, such that the file can
 * call the proper function for the selected format.
 */
template <typename format_type>
struct sequence_file_input_format_exposer : public format_type
{
public:
    // Can't use `using format_type::read_sequence_record` as it produces a hard failure in the format concept check
    // for types that do not model the format concept, i.e. don't offer the proper read_sequence_record interface.
    //!\brief Forwards to the seqan3::sequence_file_input_format::read_sequence_record interface.
    template <typename... ts>
    void read_sequence_record(ts &&... args)
    {
        format_type::read_sequence_record(std::forward<ts>(args)...);
    }
};

} // namespace seqan3::detail

namespace seqan3
{

/*!\interface seqan3::sequence_file_input_format <>
 * \brief The generic concept for sequence file in formats.
 * \ingroup io_sequence_file
 *
 * \details
 *
 * The details of this concept are only relevant to developers who wish to implement their own format.
 * The requirements for this concept are given as related functions and type traits.
 * Types that satisfy this concept are shown as "implementing this interface".
 *
 * \remark For a complete overview, take a look at \ref io_sequence_file
 */
//!\cond
template <typename t>
concept sequence_file_input_format = requires (detail::sequence_file_input_format_exposer<t> & v,
                                               std::ifstream & f,
                                               sequence_file_input_options<dna5> & options,
                                               std::streampos & position_buffer,
                                               std::vector<dna5> & seq,
                                               std::string & id,
                                               std::vector<phred42> & qual,
                                               std::vector<qualified<dna5, phred42>> & seq_qual) {
    t::file_extensions;

    { v.read_sequence_record(f, options, position_buffer, seq, id, qual) } -> std::same_as<void>;
    { v.read_sequence_record(f, options, position_buffer, seq_qual, id, seq_qual) } -> std::same_as<void>;
    {
        v.read_sequence_record(f, options, position_buffer, std::ignore, std::ignore, std::ignore)
    } -> std::same_as<void>;
};
//!\endcond

// Workaround for https://github.com/doxygen/doxygen/issues/9379
#if SEQAN3_DOXYGEN_ONLY(1) 0
template <typename t>
class sequence_file_input_format
{};
#endif

/*!\name Requirements for seqan3::sequence_file_input_format
 * \brief You can expect these **members** on all types that implement seqan3::sequence_file_input_format.
 * \memberof seqan3::sequence_file_input_format
 * \{
 */

/*!\fn void read_sequence_record(stream_type & stream,
                                 seqan3::sequence_file_input_options const & options,
                                 stream_pos_type & position_buffer,
                                 seq_type & sequence,
                                 id_type & id,
                                 qual_type & qualities)
 * \brief Read from the specified stream and back-insert into the given field buffers.
 * \tparam stream_type      Input stream, must satisfy seqan3::input_stream_over with `char`.
 * \tparam stream_pos_type  Buffer for storing the current record's file position.
 * \tparam seq_type         Type of the seqan3::field::seq input; must satisfy std::ranges::output_range
 * over a seqan3::alphabet.
 * \tparam id_type          Type of the seqan3::field::id input; must satisfy std::ranges::output_range
 * over a seqan3::alphabet.
 * \tparam qual_type        Type of the seqan3::field::qual input; must satisfy std::ranges::output_range
 * over a seqan3::writable_quality_alphabet.
 * \param[in,out] stream    The input stream to read from.
 * \param[in,out] position_buffer The buffer to store the current record's file position.
 * \param[in]     options   File specific options passed to the format.
 * \param[out]    sequence  The buffer for seqan3::field::seq input, i.e. the "sequence".
 * \param[out]    id        The buffer for seqan3::field::id input, e.g. the header line in FASTA .
 * \param[out]    qualities The buffer for seqan3::field::qual input.
 *
 * \details
 *
 * ### Additional requirements
 *
 *   * The function must also accept std::ignore as parameter for any of the fields.
 *     [This is enforced by the concept checker!]
 *   * In this case the data read for that field shall be discarded by the format.
 */
/*!\var static inline std::vector<std::string> seqan3::sequence_file_input_format::file_extensions
 * \brief The format type is required to provide a vector of all supported file extensions.
 */
//!\}

} // namespace seqan3

namespace seqan3::detail
{

/*!\brief Auxiliary value metafuncton that checks whether a type is a seqan3::type_list and all types meet
 * seqan3::sequence_file_input_format [default is false].
 * \ingroup core
 * \see seqan3::type_list_specialisationOfsequence_file_input_formats
 */
template <typename t>
constexpr bool is_type_list_of_sequence_file_input_formats_v = false;

/*!\brief Auxiliary value metafuncton that checks whether a type is a seqan3::type_list and all types meet
 * seqan3::sequence_file_input_format [overload].
 * \ingroup core
  * \see seqan3::type_list_specialisationOfsequence_file_input_formats
 */
template <typename... ts>
constexpr bool is_type_list_of_sequence_file_input_formats_v<type_list<ts...>> =
    (sequence_file_input_format<ts> && ...);

/*!\brief Auxiliary concept that checks whether a type is a seqan3::type_list and all types meet
 *        seqan3::sequence_file_input_format.
 * \ingroup core
 * \see seqan3::is_type_list_of_sequence_file_formats_v
 */
template <typename t>
concept type_list_of_sequence_file_input_formats = is_type_list_of_sequence_file_input_formats_v<t>;
} // namespace seqan3::detail
