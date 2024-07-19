// SPDX-FileCopyrightText: 2006-2024 Knut Reinert & Freie Universität Berlin
// SPDX-FileCopyrightText: 2016-2024 Knut Reinert & MPI für molekulare Genetik
// SPDX-License-Identifier: BSD-3-Clause

/*!\file
 * \brief Provides seqan3::sequence_file_output_format and auxiliary classes.
 * \author Hannes Hauswedell <hannes.hauswedell AT fu-berlin.de>
 */

#pragma once

#include <fstream>
#include <string>
#include <vector>

#include <seqan3/alphabet/nucleotide/dna5.hpp>
#include <seqan3/alphabet/quality/aliases.hpp>
#include <seqan3/alphabet/quality/phred42.hpp>
#include <seqan3/io/sequence_file/output_options.hpp>
#include <seqan3/utility/type_list/type_list.hpp>

namespace seqan3::detail
{

/*!\brief Internal class used to expose the actual format interface to write sequence records into the file.
 * \ingroup io_sequence_file
 *
 * \tparam format_type The type of the format to be exposed.
 *
 * \details
 *
 * Exposes the protected member function `write_sequence_record` from the given `format_type`, such that the file can
 * call the proper function for the selected format.
 *
 * \remark For a complete overview, take a look at \ref io_sequence_file
 */
template <typename format_type>
struct sequence_file_output_format_exposer : public format_type
{
public:
    // Can't use `using format_type::write_sequence_record` as it produces a hard failure in the format concept check
    // for types that do not model the format concept, i.e. don't offer the proper write_sequence_record interface.
    //!\brief Forwards to the seqan3::sequence_file_output_format::write_sequence_record interface.
    template <typename... ts>
    void write_sequence_record(ts &&... args)
    {
        format_type::write_sequence_record(std::forward<ts>(args)...);
    }
};

} // namespace seqan3::detail

namespace seqan3
{

/*!\interface seqan3::sequence_file_output_format <>
 * \brief The generic concept for sequence file out formats.
 * \ingroup io_sequence_file
 *
 * \details
 *
 * The details of this concept are only relevant to developers who wish to implement their own format.
 * The requirements for this concept are given as related functions and type traits.
 * Types that satisfy this concept are shown as "implementing this interface".
 */
//!\cond
template <typename t>
concept sequence_file_output_format = requires (detail::sequence_file_output_format_exposer<t> & v,
                                                std::ofstream & f,
                                                sequence_file_output_options & options,
                                                dna5_vector & seq,
                                                std::string & id,
                                                std::vector<phred42> & qual,
                                                std::vector<dna5q> & seq_qual) {
    t::file_extensions;

    { v.write_sequence_record(f, options, seq, id, qual) } -> std::same_as<void>;
    { v.write_sequence_record(f, options, std::ignore, id, std::ignore) } -> std::same_as<void>;
    // This one is required to be a valid call, but should always throw at run-time.
    { v.write_sequence_record(f, options, std::ignore, std::ignore, std::ignore) } -> std::same_as<void>;
};
//!\endcond

// Workaround for https://github.com/doxygen/doxygen/issues/9379
#if SEQAN3_DOXYGEN_ONLY(1) 0
template <typename t>
class sequence_file_output_format
{};
#endif

/*!\name Requirements for seqan3::sequence_file_output_format
 * \brief You can expect these **members** on all types that implement seqan3::sequence_file_output_format.
 * \memberof seqan3::sequence_file_output_format
 * \{
 */

/*!\fn void write_sequence_record(stream_type & stream, seqan3::sequence_file_output_options const & options, seq_type && sequence,
 *                id_type && id, qual_type && qualities)
 * \brief Write the given fields to the specified stream.
 * \tparam stream_type      Output stream, must satisfy seqan3::output_stream_over with `char`.
 * \tparam seq_type         Type of the seqan3::field::seq output; must satisfy std::ranges::output_range
 * over a seqan3::alphabet.
 * \tparam id_type          Type of the seqan3::field::id output; must satisfy std::ranges::output_range
 * over a seqan3::alphabet.
 * \tparam qual_type        Type of the seqan3::field::qual output; must satisfy std::ranges::output_range
 * over a seqan3::quality_alphabet.
 * \param[in,out] stream    The output stream to write into.
 * \param[in]     options   File specific options passed to the format.
 * \param[in]     sequence  The data for seqan3::field::seq, i.e. the "sequence".
 * \param[in]     id        The data for seqan3::field::id, e.g. the header line in FASTA.
 * \param[in]     qualities The data for seqan3::field::qual.
 *
 * \details
 *
 * ### Additional requirements
 *
 *   * The format must also accept std::ignore as parameter for any of the fields, however it shall throw an exception
 * if one of the fields required for writing the format is marked as such. [this shall be checked inside the function]
 */
/*!\var static inline std::vector<std::string> seqan3::sequence_file_output_format::file_extensions
 * \brief The format type is required to provide a vector of all supported file extensions.
 */

//!\}

} // namespace seqan3

namespace seqan3::detail
{

/*!\brief Auxiliary value metafuncton that checks whether a type is a seqan3::type_list and all types meet
 *        seqan3::sequence_file_output_format [default is false].
 * \ingroup io_sequence_file
 * \see seqan3::type_list_specialisationOfsequence_file_output_formats
 */
template <typename t>
constexpr bool is_type_list_of_sequence_file_output_formats_v = false;

/*!\brief Auxiliary value metafuncton that checks whether a type is a seqan3::type_list and all types meet
 *        seqan3::sequence_file_output_format [overload].
 * \ingroup io_sequence_file
 * \see seqan3::type_list_specialisationOfsequence_file_output_formats
 */
template <typename... ts>
constexpr bool is_type_list_of_sequence_file_output_formats_v<type_list<ts...>> =
    (sequence_file_output_format<ts> && ...);

/*!\brief Auxiliary concept that checks whether a type is a seqan3::type_list and all types meet
 *        seqan3::sequence_file_output_format.
 * \ingroup io_sequence_file
 * \see seqan3::is_type_list_of_sequence_file_formats_v
 */
template <typename t>
concept type_list_of_sequence_file_output_formats = is_type_list_of_sequence_file_output_formats_v<t>;
} // namespace seqan3::detail
