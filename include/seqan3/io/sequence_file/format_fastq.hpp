// -----------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2020, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2020, Knut Reinert & MPI für molekulare Genetik
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at: https://github.com/seqan/seqan3/blob/master/LICENSE.md
// -----------------------------------------------------------------------------------------------------

/*!\file
 * brief Provides the seqan3::format_fastq.
 * \author Hannes Hauswedell <hannes.hauswedell AT fu-berlin.de>
 */

#pragma once

#include <iterator>
#include <string>
#include <string_view>
#include <vector>

#include <seqan3/alphabet/adaptation/char.hpp>
#include <seqan3/alphabet/nucleotide/dna5.hpp>
#include <seqan3/alphabet/quality/aliases.hpp>
#include <seqan3/core/char_operations/predicate.hpp>
#include <seqan3/core/type_traits/range.hpp>
#include <seqan3/io/detail/ignore_output_iterator.hpp>
#include <seqan3/io/detail/misc.hpp>
#include <seqan3/io/sequence_file/input_format_concept.hpp>
#include <seqan3/io/sequence_file/input_options.hpp>
#include <seqan3/io/sequence_file/output_format_concept.hpp>
#include <seqan3/io/sequence_file/output_options.hpp>
#include <seqan3/io/stream/iterator.hpp>
#include <seqan3/range/detail/misc.hpp>
#include <seqan3/range/views/char_to.hpp>
#include <seqan3/range/views/istreambuf.hpp>
#include <seqan3/range/views/to_char.hpp>
#include <seqan3/range/views/take.hpp>
#include <seqan3/range/views/take_exactly.hpp>
#include <seqan3/range/views/take_line.hpp>
#include <seqan3/range/views/take_until.hpp>
#include <seqan3/std/algorithm>
#include <seqan3/std/ranges>

namespace seqan3
{

/*!\brief       The FastQ format.
 * \implements  SequenceFileFormat
 * \ingroup     sequence
 *
 * \details
 *
 * ### Introduction
 *
 * FastQ is the de-facto-standard for storing sequences together with quality information. See the
 * [article on wikipedia](https://en.wikipedia.org/wiki/FASTQ_format) for a an in-depth description of the format.
 *
 * ### fields_specialisation
 *
 * The FastQ format provides the fields seqan3::field::seq, seqan3::field::id and seqan3::field::qual; or alternatively
 * provides seqan3::field::seq_qual as a single field of sequence and quality. All three fields (or ID + SEQ_QUAL) are
 * required when writing and the sequence and qualities are required to be of the same length.
 *
 * ### Encodings
 *
 * All documented encodings for the quality string are supported (see the article above), but they are **not detected**
 * from the file. Instead, when reading the file, you have to set the respective alphabet via a traits type (see
 * seqan3::sequence_file_input_traits and the \ref quality submodule).
 *
 * ### Implementation notes
 *
 * This implementation supports the following optional features of the format:
 *
 *   * line breaks and/or other whitespace characters in any part of the sequence and/or qualities (only when reading!)
 *   * writing the ID to the `+`-line also (line is always ignored when reading)
 *
 */
class format_fastq
{
public:
    /*!\name Constructors, destructor and assignment
     * \{
     */
    format_fastq() noexcept = default; //!< Defaulted.
    format_fastq(format_fastq const &) noexcept = default; //!< Defaulted.
    format_fastq & operator=(format_fastq const &) noexcept = default; //!< Defaulted.
    format_fastq(format_fastq &&) noexcept = default; //!< Defaulted.
    format_fastq & operator=(format_fastq &&) noexcept = default; //!< Defaulted.
    ~format_fastq() noexcept = default; //!< Defaulted.
    //!\}

    //!\brief The valid file extensions for this format; note that you can modify this value.
    static inline std::vector<std::string> file_extensions
    {
        { "fastq" },
        { "fq"    }
    };

protected:
    //!\copydoc sequence_file_input_format::read_sequence_record
    template <typename stream_type,     // constraints checked by file
              typename seq_legal_alph_type, bool seq_qual_combined,
              typename seq_type,        // other constraints checked inside function
              typename id_type,
              typename qual_type>
    void read_sequence_record(stream_type                                                               & stream,
                              sequence_file_input_options<seq_legal_alph_type, seq_qual_combined> const & options,
                              seq_type                                                                  & sequence,
                              id_type                                                                   & id,
                              qual_type                                                                 & qualities)
    {
        auto stream_view = views::istreambuf(stream);
        auto stream_it = begin(stream_view);

        // cache the begin position so we write quals to the same position as seq in seq_qual case
        size_t sequence_size_before = 0;
        size_t sequence_size_after = 0;
        if constexpr (!detail::decays_to_ignore_v<seq_type>)
            sequence_size_before = size(sequence);

        /* ID */
        if (*stream_it != '@') // [[unlikely]]
        {
            throw parse_error{std::string{"Expected '@' on beginning of ID line, got: "} +
                              detail::make_printable(*stream_it)};
        }
        ++stream_it; // skip '@'

        if constexpr (!detail::decays_to_ignore_v<id_type>)
        {
            if (options.truncate_ids)
            {
                std::ranges::copy(stream_view | views::take_until_or_throw(is_cntrl || is_blank)
                                              | views::char_to<std::ranges::range_value_t<id_type>>,
                                  std::ranges::back_inserter(id));
                detail::consume(stream_view | views::take_line_or_throw);
            }
            else
            {
                std::ranges::copy(stream_view | views::take_line_or_throw
                                              | views::char_to<std::ranges::range_value_t<id_type>>,
                                  std::ranges::back_inserter(id));
            }
        }
        else
        {
            detail::consume(stream_view | views::take_line_or_throw);
        }

        /* Sequence */
        auto seq_view = stream_view | views::take_until_or_throw(is_char<'+'>)    // until 2nd ID line
                                    | std::views::filter(!is_space);           // ignore whitespace
        if constexpr (!detail::decays_to_ignore_v<seq_type>)
        {
            auto constexpr is_legal_alph = is_in_alphabet<seq_legal_alph_type>;
            std::ranges::copy(seq_view | std::views::transform([is_legal_alph] (char const c) // enforce legal alphabet
                                    {
                                        if (!is_legal_alph(c))
                                        {
                                            throw parse_error{std::string{"Encountered an unexpected letter: "} +
                                                                is_legal_alph.msg +
                                                                " evaluated to false on " +
                                                                detail::make_printable(c)};
                                        }
                                        return c;
                                    })
                                        | views::char_to<std::ranges::range_value_t<seq_type>>,         // convert to actual target alphabet
                              std::ranges::back_inserter(sequence));
            sequence_size_after = size(sequence);
        }
        else // consume, but count
        {
            auto it = begin(seq_view);
            auto it_end = end(seq_view);
            while (it != it_end)
            {
                ++it;
                ++sequence_size_after;
            }
        }

        /* 2nd ID line */
        if (*stream_it != '+') // [[unlikely]]
        {
            throw parse_error{std::string{"Expected '+' on beginning of 2nd ID line, got: "} +
                              detail::make_printable(*stream_it)};
        }
        detail::consume(stream_view | views::take_line_or_throw);

        /* Qualities */
        auto qview = stream_view | std::views::filter(!is_space)                  // this consumes trailing newline
                                 | views::take_exactly_or_throw(sequence_size_after - sequence_size_before);
        if constexpr (seq_qual_combined)
        {
            // seq_qual field implies that they are the same variable
            assert(std::addressof(sequence) == std::addressof(qualities));
            std::ranges::copy(qview | views::char_to<typename std::ranges::range_value_t<qual_type>::quality_alphabet_type>,
                              begin(qualities) + sequence_size_before);
        }
        else if constexpr (!detail::decays_to_ignore_v<qual_type>)
        {
            std::ranges::copy(qview | views::char_to<std::ranges::range_value_t<qual_type>>,
                              std::ranges::back_inserter(qualities));
        }
        else
        {
            detail::consume(qview);
        }
    }

    //!\copydoc sequence_file_output_format::write_sequence_record
    template <typename stream_type,     // constraints checked by file
              typename seq_type,        // other constraints checked inside function
              typename id_type,
              typename qual_type>
    void write_sequence_record(stream_type                     & stream,
                               sequence_file_output_options const & options,
                               seq_type                       && sequence,
                               id_type                        && id,
                               qual_type                      && qualities)
    {
        seqan3::detail::fast_ostreambuf_iterator stream_it{*stream.rdbuf()};

        // ID
        if constexpr (detail::decays_to_ignore_v<id_type>)
        {
            throw std::logic_error{"The ID field may not be set to ignore when writing FASTQ files."};
        }
        else
        {
            if (std::ranges::empty(id)) //[[unlikely]]
                throw std::runtime_error{"The ID field may not be empty when writing FASTQ files."};

            stream_it = '@';
            stream_it.write_range(id);
            stream_it.write_end_of_line(options.add_carriage_return);
        }

        // Sequence
        if constexpr (detail::decays_to_ignore_v<seq_type>)
        {
            throw std::logic_error{"The SEQ and SEQ_QUAL fields may not both be set to ignore when writing FASTQ files."};
        }
        else
        {
            if (std::ranges::empty(sequence)) //[[unlikely]]
                throw std::runtime_error{"The SEQ field may not be empty when writing FASTQ files."};

            stream_it.write_range(sequence | views::to_char);
            stream_it.write_end_of_line(options.add_carriage_return);
        }

        // 2nd ID-line
        if constexpr (!detail::decays_to_ignore_v<id_type>)
        {
            stream_it = '+';

            if (options.fastq_double_id)
                stream_it.write_range(id);

            stream_it.write_end_of_line(options.add_carriage_return);
        }

        // Quality line
        if constexpr (detail::decays_to_ignore_v<qual_type>)
        {
            throw std::logic_error{"The QUAL and SEQ_QUAL fields may not both be set to ignore when writing FASTQ files."};
        }
        else
        {
            if (std::ranges::empty(qualities)) //[[unlikely]]
                throw std::runtime_error{"The SEQ field may not be empty when writing FASTQ files."};

            if constexpr (std::ranges::sized_range<seq_type> && std::ranges::sized_range<qual_type>)
            {
                assert(std::ranges::size(sequence) == std::ranges::size(qualities));
            }

            stream_it.write_range(qualities | views::to_char);
            stream_it.write_end_of_line(options.add_carriage_return);
        }
    }
};

} // namespace seqan
