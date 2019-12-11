// -----------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2019, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2019, Knut Reinert & MPI für molekulare Genetik
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at: https://github.com/seqan/seqan3/blob/master/LICENSE.md
// -----------------------------------------------------------------------------------------------------

/*!\file
 * \brief Provides the seqan3::detail::format_sam_as_sequence_file_base which provides read and write sequence
 *        record interfaces for seqan3::format_sam and seqan3::format_bam.
 * \author René Rahn <rene.rahn AT fu-berlin.de>
 */

#pragma once

#include <optional>
#include <string>
#include <string_view>
#include <tuple>

#include <seqan3/alphabet/gap/gapped.hpp>
#include <seqan3/core/type_traits/basic.hpp>
#include <seqan3/io/alignment_file/header.hpp>
#include <seqan3/io/alignment_file/input_options.hpp>
#include <seqan3/io/alignment_file/misc.hpp>
#include <seqan3/io/alignment_file/output_options.hpp>
#include <seqan3/io/alignment_file/sam_tag_dictionary.hpp>
#include <seqan3/io/sequence_file/input_options.hpp>
#include <seqan3/io/sequence_file/output_options.hpp>
#include <seqan3/range/views/to.hpp>
#include <seqan3/std/iterator>
#include <seqan3/std/span>

namespace seqan3
{

//!\cond
// Forward declares format bam.
class format_bam;
//!\endcond

/*!\brief CRTP-base class to enable seqan3::format_sam and seqan3::format_bam as sequence file formats.
 * \ingroup alignment_file
 *
 * \tparam derived_format_t The type of the format that inherits from this CRTP-base class.
 *
 * \details
 *
 * The SAM and BAM format can also be used to store only sequence information. Accordingly, this class defines the
 * seqan3::sequence_file_input_format::read_sequence_record and
 * seqan3::sequence_file_output_format:write_sequence_record interface for the sam and bam format such that they can be
 * used inside of a seqan3::sequence_file_input and seqan3::sequence_file_output.
 */
template <typename derived_format_t>
class format_sam_as_sequence_file_base
{
private:
    //!\brief Befriend the derived format.
    friend derived_format_t;

    //!\brief Checks wether the header is required for writing sequence records into the alignment file.
    static constexpr bool header_required_for_writing = std::same_as<derived_format_t, format_bam>;

    /*!\name Constructors, destructor and assignment
     * \{
     */
    //!\brief Defaulted.
    format_sam_as_sequence_file_base() = default;
    //!\brief Defaulted.
    format_sam_as_sequence_file_base(format_sam_as_sequence_file_base const &) = default;
    //!\brief Defaulted.
    format_sam_as_sequence_file_base & operator=(format_sam_as_sequence_file_base const &) = default;
    //!\brief Defaulted.
    format_sam_as_sequence_file_base(format_sam_as_sequence_file_base &&) = default;
    //!\brief Defaulted.
    format_sam_as_sequence_file_base & operator=(format_sam_as_sequence_file_base &&) = default;
    //!\brief Defaulted.
    ~format_sam_as_sequence_file_base() = default;
    //!\}

protected:
    template <typename stream_type,     // constraints checked by file
              typename seq_legal_alph_type, bool seq_qual_combined,
              typename seq_type,        // other constraints checked inside function
              typename id_type,
              typename qual_type>
    void read_sequence_record(stream_type & stream,
                              sequence_file_input_options<seq_legal_alph_type, seq_qual_combined> const & options,
                              seq_type & sequence,
                              id_type & id,
                              qual_type & qualities);

    template <typename stream_type,     // constraints checked by file
              typename seq_type,        // other constraints checked inside function
              typename id_type,
              typename qual_type>
    void write_sequence_record(stream_type & stream,
                               sequence_file_output_options const & SEQAN3_DOXYGEN_ONLY(options),
                               seq_type && sequence,
                               id_type && id,
                               qual_type && qualities);

    //!brief Returns a reference to dummy if passed a std::ignore.
    std::string_view const & default_or(detail::ignore_t) const noexcept
    {
        return dummy;
    }

    //!brief Returns the input unchanged.
    template <typename t>
    decltype(auto) default_or(t && v) const noexcept
    {
        return std::forward<t>(v);
    }

private:
    //!\brief Converts this to the derived type.
    derived_format_t & as_derived() noexcept
    {
        return static_cast<derived_format_t &>(*this);
    }

    //!\brief Stores quality values temporarily if seq and qual information are combined (not supported by SAM yet).
    std::string tmp_qual{};

    //!\brief An empty dummy container to pass to align_format.write() such that an empty field is written.
    static constexpr std::string_view dummy{};
};

//!\copydoc seqan3::sequence_file_input_format::read_sequence_record
template <typename derived_format_t>
template <typename stream_type,     // constraints checked by file
            typename seq_legal_alph_type, bool seq_qual_combined,
            typename seq_type,        // other constraints checked inside function
            typename id_type,
            typename qual_type>
inline void format_sam_as_sequence_file_base<derived_format_t>::read_sequence_record(
                stream_type & stream,
                sequence_file_input_options<seq_legal_alph_type, seq_qual_combined> const & options,
                seq_type & sequence,
                id_type & id,
                qual_type & qualities)
{
    alignment_file_input_options<seq_legal_alph_type> align_options{};

    alignment_file_header<> default_header{};

    if constexpr (seq_qual_combined)
    {
        tmp_qual.clear();
        as_derived().read_alignment_record(
            stream, align_options, std::ignore, default_header, sequence, tmp_qual, id,
            std::ignore, std::ignore, std::ignore, std::ignore, std::ignore, std::ignore,
            std::ignore, std::ignore, std::ignore, std::ignore, std::ignore, std::ignore);

        for (auto sit = tmp_qual.begin(), dit = std::ranges::begin(sequence); sit != tmp_qual.end(); ++sit, ++dit)
            get<1>(*dit).assign_char(*sit);
    }
    else
    {
        as_derived().read_alignment_record(
            stream, align_options, std::ignore, default_header, sequence, qualities, id,
            std::ignore, std::ignore, std::ignore, std::ignore, std::ignore, std::ignore,
            std::ignore, std::ignore, std::ignore, std::ignore, std::ignore, std::ignore);
    }

    if constexpr (!detail::decays_to_ignore_v<seq_type>)
        if (std::ranges::distance(sequence) == 0)
            throw parse_error{"The sequence information must not be empty."};
    if constexpr (!detail::decays_to_ignore_v<id_type>)
        if (std::ranges::distance(id) == 0)
            throw parse_error{"The id information must not be empty."};

    if (options.truncate_ids)
        id = id | views::take_until_and_consume(is_space) | views::to<id_type>;
}

//!\copydoc seqan3::sequence_file_output_format::write_sequence_record
template <typename derived_format_t>
template <typename stream_type,     // constraints checked by file
          typename seq_type,        // other constraints checked inside function
          typename id_type,
          typename qual_type>
inline void format_sam_as_sequence_file_base<derived_format_t>::write_sequence_record(
                stream_type & stream,
                sequence_file_output_options const & SEQAN3_DOXYGEN_ONLY(options),
                seq_type && sequence,
                id_type && id,
                qual_type && qualities)
{
    using default_align_t = std::pair<std::span<gapped<char>>, std::span<gapped<char>>>;
    using default_mate_t = std::tuple<std::string_view, std::optional<int32_t>, int32_t>;
    using default_header_t = std::conditional_t<header_required_for_writing, alignment_file_header<>, detail::ignore_t>;

    alignment_file_output_options output_options{};
    default_header_t default_header{};

    as_derived().write_alignment_record(stream,
                                        output_options,
                     /*header*/         default_header,
                     /*seq*/            default_or(sequence),
                     /*qual*/           default_or(qualities),
                     /*id*/             default_or(id),
                     /*offset*/         0,
                     /*ref_seq*/        std::string_view{},
                     /*ref_id*/         std::string_view{},
                     /*ref_offset*/     -1,
                     /*align*/          default_align_t{},
                     /*cigar_vector*/   std::vector<cigar>{},
                     /*flag*/           sam_flag::none,
                     /*mapq*/           0,
                     /*mate*/           default_mate_t{},
                     /*tag_dict*/       sam_tag_dictionary{},
                     /*e_value*/        0,
                     /*bit_score*/      0);
}

} // namespace seqan3
