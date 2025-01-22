// SPDX-FileCopyrightText: 2006-2025 Knut Reinert & Freie Universität Berlin
// SPDX-FileCopyrightText: 2016-2025 Knut Reinert & MPI für molekulare Genetik
// SPDX-License-Identifier: BSD-3-Clause

/*!\file
 * \brief Provides seqan3::sam_file_output_format and auxiliary classes.
 * \author Svenja Mehringer <svenja.mehringer AT fu-berlin.de>
 */

#pragma once

#include <fstream>
#include <optional>
#include <string>
#include <vector>

#include <seqan3/alphabet/cigar/cigar.hpp>
#include <seqan3/alphabet/nucleotide/dna5.hpp>
#include <seqan3/alphabet/quality/aliases.hpp>
#include <seqan3/alphabet/quality/phred42.hpp>
#include <seqan3/io/sam_file/header.hpp>
#include <seqan3/io/sam_file/output_options.hpp>
#include <seqan3/io/sam_file/sam_flag.hpp>
#include <seqan3/io/sam_file/sam_tag_dictionary.hpp>
#include <seqan3/utility/type_list/type_list.hpp>

namespace seqan3::detail
{

/*!\brief Internal class used to expose the actual format interface to write SAM/BAM records into the file.
 * \ingroup io_sam_file
 *
 * \tparam format_type The type of the format to be exposed.
 *
 * \details
 *
 * Exposes the protected member function `write_alignment_record` and `write_header` from the given `format_type`, such
 * that the file can call the proper function for the selected format.
 *
 * \remark For a complete overview, take a look at \ref io_sam_file
 */
template <typename format_type>
struct sam_file_output_format_exposer : public format_type
{
public:
    // Can't use `using format_type::write_alignment_record` as it produces a hard failure in the format concept check
    // for types that do not model the format concept, i.e. don't offer the proper write_alignment_record interface.
    //!\brief Forwards to the seqan3::sam_file_output_format::write_alignment_record interface.
    template <typename... ts>
    void write_alignment_record(ts &&... args)
    {
        format_type::write_alignment_record(std::forward<ts>(args)...);
    }

    //!\brief Forwards to `format_type::write_header`.
    template <typename stream_t, typename header_type>
    void write_header(stream_t & stream, sam_file_output_options const & options, header_type & header)
    {
        format_type::write_header(stream, options, header);
    }
};

} // namespace seqan3::detail

namespace seqan3
{

/*!\interface seqan3::sam_file_output_format <>
 * \brief The generic concept for alignment file out formats.
 * \ingroup io_sam_file
 *
 * \details
 *
 * The details of this concept are only relevant to developers who wish to
 * implement their own format. The requirements for this concept are given as
 * related functions and type traits. Types that model this concept are shown
 * as "implementing this interface".
 */
//!\cond
template <typename t>
concept sam_file_output_format = requires (detail::sam_file_output_format_exposer<t> & v,
                                           std::ofstream & stream,
                                           sam_file_output_options & options,
                                           sam_file_header<> & header,
                                           dna5_vector & seq,
                                           std::vector<phred42> & qual,
                                           std::string & id,
                                           dna5_vector & ref_seq,
                                           std::optional<int32_t> & ref_id,
                                           std::optional<int32_t> & ref_offset,
                                           std::vector<cigar> & cigar,
                                           sam_flag & flag,
                                           uint8_t & mapq,
                                           std::tuple<std::optional<int32_t>, std::optional<int32_t>, int32_t> & mate,
                                           sam_tag_dictionary & tag_dict,
                                           double & e_value,
                                           double & bit_score) {
    t::file_extensions;

    {
        v.write_alignment_record(stream,
                                 options,
                                 header,
                                 seq,
                                 qual,
                                 id,
                                 ref_seq,
                                 ref_id,
                                 ref_offset,
                                 cigar,
                                 flag,
                                 mapq,
                                 mate,
                                 tag_dict,
                                 e_value,
                                 bit_score)
    } -> std::same_as<void>;
};
//!\endcond

// Workaround for https://github.com/doxygen/doxygen/issues/9379
#if SEQAN3_DOXYGEN_ONLY(1) 0
template <typename t>
class sam_file_output_format
{};
#endif

/*!\name Requirements for seqan3::sam_file_output_format
 * \brief You can expect these **members** on all types that implement seqan3::sam_file_output_format.
 * \memberof seqan3::sam_file_output_format
 * \{
 */

/*!\fn void write_alignment_record(stream_type                            &  stream,
                                   sam_file_output_options const          &  options,
                                   sam_file_header<>                      &  header,
                                   seq_type                               && seq,
                                   qual_type                              && qual,
                                   id_type                                && id,
                                   ref_seq_type                           && ref_seq,
                                   ref_id_type                            && ref_id,
                                   ref_offset_type                        && ref_offset,
                                   std::vector<cigar>                     &  cigar_vector,
                                   flag_type                              && flag,
                                   mapq_type                              && mapq,
                                   mate_type                              && mate,
                                   tag_dict_type                          && tag_dict,
                                   e_value_type                           && e_value,
                                   bit_score_type                         && bit_score)
 * \brief Write the given fields to the specified stream.
 * \tparam stream_type      Output stream, must model seqan3::output_stream_over with `char`.
 * \tparam seq_type         Type of the seqan3
 * \tparam id_type          Type of the seqan3
 * \tparam ref_seq_type     Type of the seqan3
 * \tparam ref_id_type      Type of the seqan3
 * \tparam ref_offset_type  Type of the seqan3
 * \tparam flag_type        Type of the seqan3
 * \tparam mapq_type        Type of the seqan3
 * \tparam qual_type        Type of the seqan3
 * \tparam mate_type        Type of the seqan3
 * \tparam tag_dict_type    Type of the seqan3
 * \tparam e_value_type     Type of the seqan3
 * \tparam bit_score_type   Type of the seqan3
 *
 * \param[in,out] stream The output stream to write into.
 * \param[in] options File specific options passed to the format.
 * \param[in] header A pointer to the header object of the file.
 * \param[in] seq The data for seqan3::field::seq, i.e. the query sequence.
 * \param[in] qual The data for seqan3::field::qual, e.g. the query quality sequence.
 * \param[in] id The data for seqan3::field::id, e.g. the read id.
 * \param[in] ref_seq The data for seqan3::field::ref_offset, i.e. the reference sequence.
 * \param[in] ref_id The data for seqan3::field::ref_id, e.g. the reference id..
 * \param[in] ref_offset The data for seqan3::field::ref_offset, i.e. the start position of the alignment in \p ref_seq.
 * \param[in] cigar_vector The data for seqan3::field::cigar, e.g. representing the alignment between query and ref.
 * \param[in] flag The data for seqan3::field::flag, e.g. the SAM mapping flag value.
 * \param[in] mapq The data for seqan3::field::mapq, e.g. the mapping quality value.
 * \param[in] mate The data for seqan3::field::mate, e.g. the mate information of paired reads.
 * \param[in] tag_dict The data for seqan3::field::tags, e.g. the optional SAM field tag dictionary.
 * \param[in] e_value The data for seqan3::field::e_value, e.g. the e-value of the alignment (BLAST).
 * \param[in] bit_score The data for seqan3::field::, e.g. the bit score of the alignment (BLAST).
 *
 */
/*!\var static inline std::vector<std::string> seqan3::sam_file_output_format::file_extensions
 * \brief The format type is required to provide a vector of all supported file extensions.
 */

//!\}

} // namespace seqan3

namespace seqan3::detail
{

/*!\brief Auxiliary value metafuncton that checks whether a type is a seqan3::type_list and all types meet
 * seqan3::sam_file_output_format [default is false].
 * \ingroup io_sam_file
 * \see seqan3::type_list_of_sam_file_output_formats
 */
template <typename t>
constexpr bool is_type_list_of_sam_file_output_formats_v = false;

/*!\brief Auxiliary value metafuncton that checks whether a type is a seqan3::type_list and all types meet
 * seqan3::sam_file_output_format [overload].
 * \ingroup io_sam_file
 * \see seqan3::type_list_of_sam_file_output_formats
 */
template <typename... ts>
constexpr bool is_type_list_of_sam_file_output_formats_v<type_list<ts...>> = (sam_file_output_format<ts> && ...);

/*!\brief Auxiliary concept that checks whether a type is a seqan3::type_list and all types meet
 * seqan3::sam_file_output_format.
 * \ingroup io_sam_file
 * \see seqan3::is_type_list_of_sam_file_output_formats_v
 */
template <typename t>
concept type_list_of_sam_file_output_formats = is_type_list_of_sam_file_output_formats_v<t>;
} // namespace seqan3::detail
