// -----------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2020, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2020, Knut Reinert & MPI für molekulare Genetik
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at: https://github.com/seqan/seqan3/blob/master/LICENSE.md
// -----------------------------------------------------------------------------------------------------

/*!\file
 * \brief Provides seqan3::alignment_file_input_format and auxiliary classes.
 * \author Svenja Mehringer <svenja.mehringer AT fu-berlin.de>
 */

#pragma once

#include <fstream>
#include <string>
#include <vector>

#include <seqan3/alphabet/cigar/cigar.hpp>
#include <seqan3/alphabet/gap/gapped.hpp>
#include <seqan3/alphabet/nucleotide/dna4.hpp>
#include <seqan3/alphabet/nucleotide/dna5.hpp>
#include <seqan3/alphabet/quality/aliases.hpp>
#include <seqan3/alphabet/quality/phred42.hpp>
#include <seqan3/core/type_list/type_list.hpp>
#include <seqan3/io/alignment_file/header.hpp>
#include <seqan3/io/alignment_file/input_options.hpp>
#include <seqan3/io/alignment_file/misc.hpp>
#include <seqan3/io/alignment_file/sam_tag_dictionary.hpp>

namespace seqan3::detail
{

/*!\brief Internal class used to expose the actual format interface to read alignment records from the file.
 * \ingroup alignment_file
 *
 * \tparam format_type The type of the format to be exposed.
 *
 * \details
 *
 * Exposes the protected member function `read_alignment_record` from the given `format_type`, such that the file can
 * call the proper function for the selected format.
 */
template <typename format_type>
struct alignment_file_input_format_exposer : public format_type
{
public:
    // Can't use `using format_type::read_alignment_record` as it produces a hard failure in the format concept check
    // for types that do not model the format concept, i.e. don't offer the proper read_alignment_record interface.
    //!\brief Forwards to the seqan3::alignment_file_input_format::read_alignment_record interface.
    template <typename ...ts>
    void read_alignment_record(ts && ...args)
    {
        format_type::read_alignment_record(std::forward<ts>(args)...);
    }
};

} // namespace seqan3::detail

namespace seqan3
{

/*!\interface seqan3::alignment_file_input_format <>
 * \brief The generic concept for alignment file input formats.
 * \ingroup alignment_file
 *
 * \details
 *
 * The details of this concept are only relevant to developers who wish to implement their own format.
 * The requirements for this concept are given as related functions and type traits.
 * Types that satisfy this concept are shown as "implementing this interface".
 */
//!\cond
template <typename t>
SEQAN3_CONCEPT alignment_file_input_format =
    requires (detail::alignment_file_input_format_exposer<t>                      & v,
              std::ifstream                                                       & stream,
              alignment_file_input_options<dna5>                                  & options,
              std::vector<dna5_vector>                                            & ref_sequences,
              alignment_file_header<>                                             & header,
              dna5_vector                                                         & seq,
              std::vector<phred42>                                                & qual,
              std::string                                                         & id,
              int32_t                                                             & offset,
              dna5_vector                                                         & ref_seq,
              std::optional<int32_t>                                              & ref_id,
              std::optional<int32_t>                                              & ref_offset,
              std::pair<std::vector<gapped<dna4>>, std::vector<gapped<dna4>>>     & align,
              std::vector<cigar>                                                  & cigar,
              sam_flag                                                            & flag,
              uint8_t                                                             & mapq,
              std::tuple<std::optional<int32_t>, std::optional<int32_t>, int32_t> & mate,
              sam_tag_dictionary                                                  & tag_dict,
              double                                                              & e_value,
              double                                                              & bit_score)
{
    t::file_extensions;
    // std::same_as<decltype(t::file_extensions), std::vector<std::string>>;

    { v.read_alignment_record(stream,
                              options,
                              ref_sequences,
                              header,
                              seq,
                              qual,
                              id,
                              offset,
                              ref_seq,
                              ref_id,
                              ref_offset,
                              align,
                              cigar,
                              flag,
                              mapq,
                              mate,
                              tag_dict,
                              e_value,
                              bit_score)};

    { v.read_alignment_record(stream,
                              options,
                              std::ignore,
                              header,
                              std::ignore,
                              std::ignore,
                              std::ignore,
                              std::ignore,
                              std::ignore,
                              std::ignore,
                              std::ignore,
                              std::ignore,
                              std::ignore,
                              std::ignore,
                              std::ignore,
                              std::ignore,
                              std::ignore,
                              std::ignore,
                              std::ignore)};
};
//!\endcond

/*!\name Requirements for seqan3::alignment_file_input_format
 * \brief You can expect these **members** on all types that implement seqan3::alignment_file_input_format.
 * \memberof seqan3::alignment_file_input_format
 * \{
 */

/*!\fn void read_alignment_record(stream_type & stream,
 *                                alignment_file_input_options<seq_legal_alph_type> const & options,
 *                                ref_seqs_type & ref_seqs,
 *                                header_type & header,
 *                                seq_type & seq,
 *                                qual_type & qual,
 *                                id_type & id,
 *                                offset_type & offset,
 *                                ref_seq_type & ref_seq,
 *                                ref_id_type & ref_id,
 *                                ref_offset_type & ref_offset,
 *                                align_type & align,
 *                                cigar_type & cigar_vector,
 *                                flag_type & flag,
 *                                mapq_type & mapq,
 *                                mate_type & mate,
 *                                tag_dict_type & tag_dict,
 *                                e_value_type & e_value,
 *                                bit_score_type & bit_score)
 * \brief Read from the specified stream and back-insert into the given field buffers.
 * \tparam stream_type        The input stream type; Must be derived from std::ostream.
 * \tparam ref_seqs_type      e.g. std::deque<ref_sequence_type> or decltype(std::ignore).
 * \tparam seq_type           Type of the seqan3::field::seq input (see seqan3::alignment_file_input_traits).
 * \tparam qual_type          Type of the seqan3::field::qual input (see seqan3::alignment_file_input_traits).
 * \tparam id_type            Type of the seqan3::field::id input (see seqan3::alignment_file_input_traits).
 * \tparam offset_type        Type of the seqan3::field::offset input (see seqan3::alignment_file_input_traits).
 * \tparam ref_seq_type       Type of the seqan3::field::ref_seq input (see seqan3::alignment_file_input_traits).
 * \tparam ref_id_type        Type of the seqan3::field::ref_id input (see seqan3::alignment_file_input_traits).
 * \tparam ref_offset_type    Type of the seqan3::field::ref_offset input (see seqan3::alignment_file_input_traits).
 * \tparam align_type         Type of the seqan3::field::alignment input (see seqan3::alignment_file_input_traits).
 * \tparam cigar_type         Type of the seqan3::field::cigar input (a std::vector<cigar> or std::ignore).
 * \tparam flag_type          Type of the seqan3::field::flag input (see seqan3::alignment_file_input_traits).
 * \tparam mapq_type          Type of the seqan3::field::mapq input (see seqan3::alignment_file_input_traits).
 * \tparam mate_type          std::tuple<ref_id_type, ref_offset_type, int32_t> or decltype(std::ignore).
 * \tparam tag_dict_type      seqan3::sam_tag_dictionary or decltype(std::ignore).
 * \tparam e_value_type       Type of the seqan3::field::evalue input (see seqan3::alignment_file_input_traits).
 * \tparam bit_score_type     Type of the seqan3::field::bit_score input (see seqan3::alignment_file_input_traits).
 *
 * \param[in,out] stream      The input stream to read from.
 * \param[in]     options     File specific options passed to the format.
 * \param[out]    ref_seqs    The reference sequences to the corresponding alignments.
 * \param[out]    header      A pointer to the seqan3::alignment_file_header object.
 * \param[out]    seq         The buffer for seqan3::field::seq input.
 * \param[out]    qual        The buffer for seqan3::field::qual input.
 * \param[out]    id          The buffer for seqan3::field::id input.
 * \param[out]    offset      The buffer for seqan3::field::offset input.
 * \param[out]    ref_seq     The buffer for seqan3::field::ref_seq input.
 * \param[out]    ref_id      The buffer for seqan3::field::ref_id input.
 * \param[out]    ref_offset  The buffer for seqan3::field::ref_offset input.
 * \param[out]    align       The buffer for seqan3::field::alignment input.
 * \param[out]    cigar_vector The buffer for seqan3::field::cigar input.
 * \param[out]    flag        The buffer for seqan3::field::flag input.
 * \param[out]    mapq        The buffer for seqan3::field::mapq input.
 * \param[out]    mate        The buffer for seqan3::field::mate input.
 * \param[out]    tag_dict    The buffer for seqan3::field::tags input.
 * \param[out]    e_value     The buffer for seqan3::field::evalue input.
 * \param[out]    bit_score   The buffer for seqan3::field::bit_score input.
 *
 * \details
 *
 * ### Additional requirements
 *
 *   * The function must also accept std::ignore as parameter for any of the fields,
 *     except stream, options and header. [This is enforced by the concept checker!]
 *   * In this case the data read for that field shall be discarded by the format.
 */
 /*!\var static inline std::vector<std::string> seqan3::alignment_file_input_format::file_extensions
 * \brief The format type is required to provide a vector of all supported file extensions.
 */
//!\}

} // namespace seqan3

namespace seqan3::detail
{

/*!\brief Auxiliary value metafuncton that checks whether a type is a seqan3::type_list and all types meet
 * seqan3::alignment_file_input_format [default is false].
 * \ingroup core
 * \see seqan3::type_list_specialisationOfalignment_file_input_formats
 */
template <typename t>
constexpr bool is_type_list_of_alignment_file_input_formats_v = false;

/*!\brief Auxiliary value metafuncton that checks whether a type is a seqan3::type_list and all types meet
 * seqan3::alignment_file_input_format [overload].
 * \ingroup core
  * \see seqan3::type_list_specialisationOfalignment_file_input_formats
 */
template <typename ...ts>
constexpr bool is_type_list_of_alignment_file_input_formats_v<type_list<ts...>> =
    (alignment_file_input_format<ts> && ...);

/*!\brief Auxiliary concept that checks whether a type is a seqan3::type_list and all types meet
 *        seqan3::alignment_file_input_format.
 * \ingroup core
 * \see seqan3::is_type_list_of_alignment_file_formats_v
 */
template <typename t>
SEQAN3_CONCEPT type_list_of_alignment_file_input_formats = is_type_list_of_alignment_file_input_formats_v<t>;

} // namespace seqan3::detail
