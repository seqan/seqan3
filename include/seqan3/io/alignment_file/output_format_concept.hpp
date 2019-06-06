// -----------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2019, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2019, Knut Reinert & MPI für molekulare Genetik
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at: https://github.com/seqan/seqan3/blob/master/LICENSE.md
// -----------------------------------------------------------------------------------------------------

/*!\file
 * \brief Provides seqan3::AlignmentFileOutputFormat and auxiliary classes.
 * \author Svenja Mehringer <svenja.mehringer AT fu-berlin.de>
 */

#pragma once

#include <fstream>
#include <string>
#include <vector>

#include <seqan3/alphabet/gap/gapped.hpp>
#include <seqan3/alphabet/nucleotide/dna5.hpp>
#include <seqan3/alphabet/quality/aliases.hpp>
#include <seqan3/alphabet/quality/phred42.hpp>
#include <seqan3/core/type_list.hpp>
#include <seqan3/io/alignment_file/header.hpp>
#include <seqan3/io/alignment_file/output_options.hpp>
#include <seqan3/io/alignment_file/sam_tag_dictionary.hpp>

namespace seqan3::detail
{

//!\brief The alignment file output format base class.
template <typename t>
class alignment_file_output_format
{};

} // namespace seqan3::detail

namespace seqan3
{

/*!\interface seqan3::AlignmentFileOutputFormat <>
 * \brief The generic concept for alignment file out formats.
 * \ingroup alignment_file
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
SEQAN3_CONCEPT AlignmentFileOutputFormat =
    requires (detail::alignment_file_output_format<t>                              & v,
              std::ofstream                                                        & stream,
              alignment_file_output_options                                        & options,
              alignment_file_header<>                                              & header,
              dna5_vector                                                          & seq,
              std::vector<phred42>                                                 & qual,
              std::string                                                          & id,
              int32_t                                                              & offset,
              dna5_vector                                                          & ref_seq,
              std::optional<int32_t>                                               & ref_id,
              std::optional<int32_t>                                               & ref_offset,
              std::pair<std::vector<gapped<dna4>>, std::vector<gapped<dna4>>>      & align,
              uint16_t                                                             & flag,
              uint8_t                                                              & mapq,
              std::tuple<std::optional<int32_t>, std::optional<int32_t>, int32_t>  & mate,
              sam_tag_dictionary                                                   & tag_dict,
              double                                                               & e_value,
              double                                                               & bit_score)
{
    t::file_extensions;

    { v.write(stream,
              options,
              header,
              seq,
              qual,
              id,
              offset,
              ref_seq,
              ref_id,
              ref_offset,
              align,
              flag,
              mapq,
              mate,
              tag_dict,
              e_value,
              bit_score
              ) } -> void;
};
//!\endcond

/*!\name Requirements for seqan3::AlignmentFileOutputFormat
 * \brief You can expect these **members** on all types that implement seqan3::AlignmentFileOutputFormat.
 * \memberof seqan3::AlignmentFileOutputFormat
 * \{
 */

/*!\fn void write(stream_type                            &  stream,
                  alignment_file_output_options const    &  options,
                  alignment_file_header<>                & header,
                  seq_type                               && seq,
                  qual_type                              && qual,
                  id_type                                && id,
                  offset_type                            && offset,
                  ref_seq_type                           && ref_seq,
                  ref_id_type                            && ref_id,
                  ref_offset_type                        && ref_offset,
                  align_type                             && align,
                  flag_type                              && flag,
                  mapq_type                              && mapq,
                  mate_type                              && mate,
                  tag_dict_type                          && tag_dict,
                  e_value_type                           && e_value,
                  bit_score_type                         && bit_score)
 * \brief Write the given fields to the specified stream.
 * \tparam stream_type      Output stream, must model seqan3::OStream with `char`.
 * \tparam seq_type         Type of the seqan3
 * \tparam id_type          Type of the seqan3
 * \tparam offset_type      Type of the seqan3
 * \tparam ref_seq_type     Type of the seqan3
 * \tparam ref_id_type      Type of the seqan3
 * \tparam ref_offset_type  Type of the seqan3
 * \tparam align_type       Type of the seqan3
 * \tparam flag_type        Type of the seqan3
 * \tparam mapq_type        Type of the seqan3
 * \tparam qual_type        Type of the seqan3
 * \tparam mate_type        Type of the seqan3
 * \tparam tag_dict_type    Type of the seqan3
 * \tparam e_value_type     Type of the seqan3
 * \tparam bit_score_type   Type of the seqan3
 *
 * \param[in,out] stream     The output stream to write into.
 * \param[in]     options    File specific options passed to the format.
 * \param[in]     header     A pointer to the header object of the file.
 * \param[in]     seq        The data for seqan3::field::SEQ, i.e. the query sequence.
 * \param[in]     qual       The data for seqan3::field::QUAL, e.g. the query quality sequence.
 * \param[in]     id         The data for seqan3::field::ID, e.g. the read id.
 * \param[in]     offset     The data for seqan3::field::OFFSET, i.e. the start position of the alignment in \p seq.
 * \param[in]     ref_seq    The data for seqan3::field::REF_OFFSET, i.e. the reference sequence.
 * \param[in]     ref_id     The data for seqan3::field::REF_ID, e.g. the reference id..
 * \param[in]     ref_offset The data for seqan3::field::REF_OFFSET, i.e. the start position of the alignment in \p ref_seq.
 * \param[in]     align      The data for seqan3::field::ALIGN, e.g. the alignment between query and ref.
 * \param[in]     flag       The data for seqan3::field::FLAG, e.g. the SAM mapping flag value.
 * \param[in]     mapq       The data for seqan3::field::MAPQ, e.g. the mapping quality value.
 * \param[in]     mate       The data for seqan3::field::MATE, e.g. the mate information of paired reads.
 * \param[in]     tag_dict   The data for seqan3::field::TAGS, e.g. the optional SAM field tag dictionary.
 * \param[in]     e_value    The data for seqan3::field::E_VALUE, e.g. the e-value of the alignment (BLAST).
 * \param[in]     bit_score  The data for seqan3::field::, e.g. the bit score of the alignment (BLAST).
 *
 */
/*!\var static inline std::vector<std::string> seqan3::AlignmentFileOutputFormat::file_extensions
 * \brief The format type is required to provide a vector of all supported file extensions.
 */

//!\}

} // namespace seqan3

namespace seqan3::detail
{

/*!\brief Auxiliary value metafuncton that checks whether a type is a seqan3::type_list and all types meet
 * seqan3::AlignmentFileOutputFormat [default is false].
 * \ingroup core
 * \see seqan3::TypeListOfAlignmentFileOutputFormats
 */
template <typename t>
constexpr bool is_type_list_of_alignment_file_output_formats_v = false;

/*!\brief Auxiliary value metafuncton that checks whether a type is a seqan3::type_list and all types meet
 * seqan3::AlignmentFileOutputFormat [overload].
 * \ingroup core
 * \see seqan3::TypeListOfAlignmentFileOutputFormats
 */
template <typename ... ts>
constexpr bool is_type_list_of_alignment_file_output_formats_v<type_list<ts...>>
                = (AlignmentFileOutputFormat<ts> && ...);

/*!\brief Auxiliary concept that checks whether a type is a seqan3::type_list and all types meet
 * seqan3::AlignmentFileOutputFormat.
 * \ingroup core
 * \see seqan3::is_type_list_of_alignment_file_formats_v
 */
template <typename t>
SEQAN3_CONCEPT TypeListOfAlignmentFileOutputFormats = is_type_list_of_alignment_file_output_formats_v<t>;
} // namespace seqan3::detail
