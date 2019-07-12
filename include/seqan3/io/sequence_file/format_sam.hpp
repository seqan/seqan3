// -----------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2019, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2019, Knut Reinert & MPI für molekulare Genetik
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at: https://github.com/seqan/seqan3/blob/master/LICENSE.md
// -----------------------------------------------------------------------------------------------------

/*!\file
 * \brief Provides the seqan3::format_sam tag and the seqan3::sequence_file_input_format and
 *        seqan3::sequence_file_output_format specialisation for this tag.
 * \author Mitra Darvish <mitra.darvish AT fu-berlin.de>
 */

#pragma once

#include <iterator>
#include <string>
#include <string_view>
#include <vector>

#include <seqan3/io/alignment_file/format_sam.hpp>
#include <seqan3/io/sequence_file/input_format_concept.hpp>
#include <seqan3/io/sequence_file/input_options.hpp>
#include <seqan3/io/sequence_file/output_format_concept.hpp>
#include <seqan3/io/sequence_file/output_options.hpp>

namespace seqan3::detail
{

/*!\brief The seqan3::sequence_file_input_format specialisation that handles formatted SAM input.
 * \ingroup sequence
 *
 * \details
 *
 * ### Introduction
 *
 * The SAM format is commonly used to store pairwise alignment information between a query sequence and its
 * reference sequence, e.g. a read mapping result. Some people also use the SAM format as plain storage for
 * sequences (and qualities) and in some cases the original sequence files are no longer available. The
 * seqan3::format_sam allows using SAM files in this manner and provides easy convertibility
 * from/to FASTQ; but there is no access to the alignment information stored in SAM files. Use
 * seqan3::format_sam if you are interested in the alignment.
 *
 * See seqan3::format_sam
 *
 */
template <>
class sequence_file_input_format<format_sam>
{
public:
    //!\brief Exposes the format tag that this class is specialised with.
    using format_tag = format_sam;

    /*!\name Constructors, destructor and assignment
     * \{
     */
    sequence_file_input_format()                                                = default; //!< Defaulted.
    //!\brief Copy construction is explicitly deleted, because you can't have multiple access to the same file.
    sequence_file_input_format(sequence_file_input_format const &)              = delete;
    //!\brief Copy assignment is explicitly deleted, because you can't have multiple access to the same file.
    sequence_file_input_format & operator=(sequence_file_input_format const &)  = delete;
    sequence_file_input_format(sequence_file_input_format &&)                   = default; //!< Defaulted.
    sequence_file_input_format & operator=(sequence_file_input_format &&)       = default; //!< Defaulted.
    ~sequence_file_input_format()                                               = default; //!< Defaulted.
    //!\}

    //!\copydoc SequenceFileInputFormat::read
    template <typename stream_type,     // constraints checked by file
              typename seq_legal_alph_type, bool seq_qual_combined,
              typename seq_type,        // other constraints checked inside function
              typename id_type,
              typename qual_type>
    void read(stream_type                                                               & stream,
              sequence_file_input_options<seq_legal_alph_type, seq_qual_combined> const & SEQAN3_DOXYGEN_ONLY(options),
              seq_type                                                                  & sequence,
              id_type                                                                   & id,
              qual_type                                                                 & qualities)
    {
        alignment_file_input_options<seq_legal_alph_type> align_options;

        if constexpr (seq_qual_combined)
        {
            tmp_qual.clear();
            align_format.read(stream, align_options, std::ignore, default_header, sequence, tmp_qual, id,
                              std::ignore, std::ignore, std::ignore, std::ignore, std::ignore,
                              std::ignore, std::ignore, std::ignore, std::ignore, std::ignore, std::ignore);

            for (auto sit = tmp_qual.begin(), dit = std::ranges::begin(sequence); sit != tmp_qual.end(); ++sit, ++dit)
                get<1>(*dit).assign_char(*sit);
        }
        else
        {
            align_format.read(stream, align_options, std::ignore, default_header, sequence, qualities, id,
                              std::ignore, std::ignore, std::ignore, std::ignore, std::ignore,
                              std::ignore, std::ignore, std::ignore, std::ignore, std::ignore, std::ignore);
        }

        if constexpr (!detail::decays_to_ignore_v<seq_type>)
            if (std::distance(std::ranges::begin(sequence), std::ranges::end(sequence)) == 0)
                throw format_error{"The sequence information must not be empty."};
        if constexpr (!detail::decays_to_ignore_v<id_type>)
            if (std::distance(std::ranges::begin(id), std::ranges::end(id)) == 0)
                throw format_error{"The sequence information must not be empty."};
    }

private:
    //!\brief An instance of the alignment format to read formatted SAM input.
    alignment_file_input_format<format_sam> align_format{};

    //!\brief The default header for the alignment format.
    alignment_file_header<> default_header{};

    //!\brief Stores quality values temporarily if seq and qual information are combined (not supported by SAM yet).
    std::string tmp_qual{};
};

//!\brief The seqan3::sequence_file_output_format specialisation that can write formatted SAM.
//!\ingroup sequence
template <>
class sequence_file_output_format<format_sam>
{
public:
    //!\brief Exposes the format tag that this class is specialised with.
    using format_tag = format_sam;

    /*!\name Constructors, destructor and assignment
     * \{
     */
    sequence_file_output_format()                                                noexcept = default; //!< Defaulted.
    //!\brief Copy construction is explicitly deleted, because you can't have multiple access to the same file.
    sequence_file_output_format(sequence_file_output_format const &)                      = delete;
    //!\brief Copy assignment is explicitly deleted, because you can't have multiple access to the same file.
    sequence_file_output_format & operator=(sequence_file_output_format const &)          = delete;
    sequence_file_output_format(sequence_file_output_format &&)                  noexcept = default; //!< Defaulted.
    sequence_file_output_format & operator=(sequence_file_output_format &&)      noexcept = default; //!< Defaulted.
    ~sequence_file_output_format()                                               noexcept = default; //!< Defaulted.
    //!\}

    //!\copydoc SequenceFileOutputFormat::write
    template <typename stream_type,     // constraints checked by file
              typename seq_type,        // other constraints checked inside function
              typename id_type,
              typename qual_type>
    void write(stream_type                        & stream,
               sequence_file_output_options const & SEQAN3_DOXYGEN_ONLY(options),
               seq_type                           && sequence,
               id_type                            && id,
               qual_type                          && qualities)
    {
        using default_align_t = std::pair<std::span<gapped<char>>, std::span<gapped<char>>>;
        using default_mate_t  = std::tuple<std::string_view, std::optional<int32_t>, int32_t>;

        alignment_file_output_options output_options;

        align_format.write(stream, output_options, std::ignore,
                           default_or(sequence), default_or(qualities), default_or(id),
                           0, std::string_view{}, std::string_view{}, -1, default_align_t{}, 0, 0,
                           default_mate_t{}, sam_tag_dictionary{}, 0, 0);
    }

private:
    //!\brief An instance of the alignment format to read formatted SAM input.
    alignment_file_output_format<format_sam> align_format{};

    //!\brief An empty dummy container to pass to align_format.write() such that an empty field is written.
    static constexpr std::string_view dummy{};

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
};

} // namespace seqan3::detail
