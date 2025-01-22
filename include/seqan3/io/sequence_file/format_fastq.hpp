// SPDX-FileCopyrightText: 2006-2025 Knut Reinert & Freie Universität Berlin
// SPDX-FileCopyrightText: 2016-2025 Knut Reinert & MPI für molekulare Genetik
// SPDX-License-Identifier: BSD-3-Clause

/*!\file
 * brief Provides the seqan3::format_fastq.
 * \author Hannes Hauswedell <hannes.hauswedell AT fu-berlin.de>
 */

#pragma once

#include <algorithm>
#include <iterator>
#include <ranges>
#include <string>
#include <string_view>
#include <vector>

#include <seqan3/alphabet/adaptation/char.hpp>
#include <seqan3/alphabet/concept.hpp>
#include <seqan3/alphabet/nucleotide/dna5.hpp>
#include <seqan3/alphabet/quality/aliases.hpp>
#include <seqan3/alphabet/views/char_to.hpp>
#include <seqan3/alphabet/views/to_char.hpp>
#include <seqan3/core/range/detail/misc.hpp>
#include <seqan3/core/range/type_traits.hpp>
#include <seqan3/io/detail/ignore_output_iterator.hpp>
#include <seqan3/io/detail/misc.hpp>
#include <seqan3/io/sequence_file/input_format_concept.hpp>
#include <seqan3/io/sequence_file/input_options.hpp>
#include <seqan3/io/sequence_file/output_format_concept.hpp>
#include <seqan3/io/sequence_file/output_options.hpp>
#include <seqan3/io/stream/detail/fast_ostreambuf_iterator.hpp>
#include <seqan3/io/views/detail/istreambuf_view.hpp>
#include <seqan3/io/views/detail/take_exactly_view.hpp>
#include <seqan3/io/views/detail/take_line_view.hpp>
#include <seqan3/io/views/detail/take_until_view.hpp>
#include <seqan3/utility/char_operations/predicate.hpp>
#include <seqan3/utility/concept.hpp>
#include <seqan3/utility/detail/type_name_as_string.hpp>

namespace seqan3
{

/*!\brief The FASTQ format.
 * \implements SequenceFileFormat
 * \ingroup io_sequence_file
 *
 * \details
 *
 * ### Introduction
 *
 * FASTQ is the de-facto-standard for storing sequences together with quality information. See the
 * [article on wikipedia](https://en.wikipedia.org/wiki/FASTQ_format) for a an in-depth description of the format.
 *
 * ### fields_specialisation
 *
 * The FASTQ format provides the fields seqan3::field::seq, seqan3::field::id and seqan3::field::qual. All three fields
 * are required when writing and the sequence and qualities are required to be of the same length.
 *
 * ### Encodings
 *
 * All documented encodings for the quality string are supported (see the article above), but they are **not detected**
 * from the file. Instead, when reading the file, you have to set the respective alphabet via a traits type (see
 * seqan3::sequence_file_input_traits and the \ref alphabet_quality submodule).
 *
 * ### Implementation notes
 *
 * This implementation supports the following optional features of the format:
 *
 *   * line breaks and/or other whitespace characters in any part of the sequence and/or qualities (only when reading!)
 *   * writing the ID to the `+`-line also (line is always ignored when reading)
 *
 * \remark For a complete overview, take a look at \ref io_sequence_file
 */
class format_fastq
{
public:
    /*!\name Constructors, destructor and assignment
     * \{
     */
    format_fastq() noexcept = default;                                 //!< Defaulted.
    format_fastq(format_fastq const &) noexcept = default;             //!< Defaulted.
    format_fastq & operator=(format_fastq const &) noexcept = default; //!< Defaulted.
    format_fastq(format_fastq &&) noexcept = default;                  //!< Defaulted.
    format_fastq & operator=(format_fastq &&) noexcept = default;      //!< Defaulted.
    ~format_fastq() noexcept = default;                                //!< Defaulted.

    //!\}

    //!\brief The valid file extensions for this format; note that you can modify this value.
    static inline std::vector<std::string> file_extensions{{"fastq"}, {"fq"}};

protected:
    //!\copydoc sequence_file_input_format::read_sequence_record
    template <typename stream_type, // constraints checked by file
              typename seq_legal_alph_type,
              typename stream_pos_type,
              typename seq_type, // other constraints checked inside function
              typename id_type,
              typename qual_type>
    void read_sequence_record(stream_type & stream,
                              sequence_file_input_options<seq_legal_alph_type> const & options,
                              stream_pos_type & position_buffer,
                              seq_type & sequence,
                              id_type & id,
                              qual_type & qualities)
    {
        // Store current position in buffer
        // Must happen before constructing the view.
        // With libc++, tellg invalidates the I/O buffer.
        position_buffer = stream.tellg();

        auto stream_view = detail::istreambuf(stream);
        auto stream_it = std::ranges::begin(stream_view);

        // cache the begin position so we write quals to the same position as seq in seq_qual case
        size_t sequence_size_before = 0;
        size_t sequence_size_after = 0;
        if constexpr (!detail::decays_to_ignore_v<seq_type>)
            sequence_size_before = size(sequence);

        /* ID */
        if (*stream_it != '@') // [[unlikely]]
        {
            throw parse_error{std::string{"Expected '@' on beginning of ID line, got: "}
                              + detail::make_printable(*stream_it)};
        }
        ++stream_it; // skip '@'

#if SEQAN3_WORKAROUND_VIEW_PERFORMANCE // can't have nice things :'(
        auto e = std::ranges::end(stream_view);
        if constexpr (!detail::decays_to_ignore_v<id_type>)
        {
            if (options.truncate_ids)
            {
                for (; (stream_it != e) && (!(is_cntrl || is_blank))(*stream_it); ++stream_it)
                {
                    if constexpr (builtin_character<std::ranges::range_value_t<id_type>>)
                        id.push_back(*stream_it);
                    else
                        id.push_back(assign_char_to(*stream_it, std::ranges::range_value_t<id_type>{}));
                }
                for (; (stream_it != e) && (!is_char<'\n'>)(*stream_it); ++stream_it)
                {}
            }
            else
            {
                for (; (stream_it != e) && (!is_char<'\n'>)(*stream_it); ++stream_it)
                {
                    if constexpr (builtin_character<std::ranges::range_value_t<id_type>>)
                        id.push_back(*stream_it);
                    else
                        id.push_back(assign_char_to(*stream_it, std::ranges::range_value_t<id_type>{}));
                }
            }
        }
        else
        {
            for (; (stream_it != e) && (!is_char<'\n'>)(*stream_it); ++stream_it)
            {}
        }

        if (stream_it == e)
        {
            throw unexpected_end_of_input{"Expected end of ID-line, got end-of-file."};
        }
        ++stream_it; // skip newline

        /* Sequence */
        if constexpr (!detail::decays_to_ignore_v<seq_type>)
        {
            for (; (stream_it != e) && (!is_char<'+'>)(*stream_it); ++stream_it)
            {
                if ((!is_space)(*stream_it))
                {
                    if constexpr (builtin_character<std::ranges::range_value_t<seq_type>>)
                    {
                        sequence.push_back(*stream_it);
                    }
                    else
                    {
                        if (!char_is_valid_for<seq_legal_alph_type>(*stream_it))
                        {
                            throw parse_error{std::string{"Encountered bad letter for seq: "}
                                              + detail::make_printable(*stream_it)};
                        }
                        sequence.push_back(assign_char_to(*stream_it, std::ranges::range_value_t<seq_type>{}));
                    }
                }
            }
            sequence_size_after = size(sequence);
        }
        else // consume, but count
        {
            for (; (stream_it != e) && (!is_char<'+'>)(*stream_it); ++stream_it)
                if ((!is_space)(*stream_it))
                    ++sequence_size_after;
        }

        /* 2nd ID line */
        if (stream_it == e)
            throw unexpected_end_of_input{"Expected second ID-line, got end-of-file."};

        if (*stream_it != '+')
        {
            throw parse_error{std::string{"Expected '+' on beginning of 2nd ID line, got: "}
                              + detail::make_printable(*stream_it)};
        }

        for (; (stream_it != e) && (!is_char<'\n'>)(*stream_it); ++stream_it)
        {}

        if (stream_it == e)
            throw unexpected_end_of_input{"Expected end of second ID-line, got end-of-file."};

        ++stream_it;

        /* Qualities */
        if constexpr (!detail::decays_to_ignore_v<qual_type>)
        {
            while (sequence_size_after > sequence_size_before)
            {
                if (stream_it == e)
                    throw unexpected_end_of_input{"Expected qualities, got end-of-file."};

                if ((!is_space)(*stream_it))
                {
                    --sequence_size_after;
                    if constexpr (builtin_character<std::ranges::range_value_t<qual_type>>)
                    {
                        qualities.push_back(*stream_it);
                    }
                    else
                    {
                        if (!char_is_valid_for<std::ranges::range_value_t<qual_type>>(*stream_it))
                        {
                            throw parse_error{std::string{"Encountered bad letter for qual: "}
                                              + detail::make_printable(*stream_it)};
                        }
                        qualities.push_back(assign_char_to(*stream_it, std::ranges::range_value_t<qual_type>{}));
                    }
                }
                ++stream_it;
            }
        }
        else // consume
        {
            while (sequence_size_after > sequence_size_before)
            {
                if (stream_it == e)
                    throw unexpected_end_of_input{"File ended before expected number of qualities could be read."};

                if ((!is_space)(*stream_it))
                    --sequence_size_after;
                ++stream_it;
            }
        }

        if (stream_it != e)
        {
            if ((!is_char<'\n'>)(*stream_it))
                throw parse_error{"Qualitites longer than sequence."};
            else
                ++stream_it;
        }

#else // ↑↑↑ WORKAROUND | ORIGINAL ↓↓↓

        if constexpr (!detail::decays_to_ignore_v<id_type>)
        {
            if (options.truncate_ids)
            {
                std::ranges::copy(stream_view | detail::take_until_or_throw(is_cntrl || is_blank)
                                      | views::char_to<std::ranges::range_value_t<id_type>>,
                                  std::back_inserter(id));
                detail::consume(stream_view | detail::take_line_or_throw);
            }
            else
            {
                std::ranges::copy(stream_view | detail::take_line_or_throw
                                      | views::char_to<std::ranges::range_value_t<id_type>>,
                                  std::back_inserter(id));
            }
        }
        else
        {
            detail::consume(stream_view | detail::take_line_or_throw);
        }

        /* Sequence */
        auto seq_view = stream_view | detail::take_until_or_throw(is_char<'+'>) // until 2nd ID line
                      | std::views::filter(!is_space);                          // ignore whitespace
        if constexpr (!detail::decays_to_ignore_v<seq_type>)
        {
            constexpr auto is_legal_alph = char_is_valid_for<seq_legal_alph_type>;
            std::ranges::copy(
                seq_view
                    | std::views::transform(
                        [is_legal_alph](char const c) // enforce legal alphabet
                        {
                            if (!is_legal_alph(c))
                            {
                                throw parse_error{std::string{"Encountered an unexpected letter: "}
                                                  + "char_is_valid_for<"
                                                  + detail::type_name_as_string<seq_legal_alph_type>
                                                  + "> evaluated to false on " + detail::make_printable(c)};
                            }
                            return c;
                        })
                    | views::char_to<std::ranges::range_value_t<seq_type>>, // convert to actual target alphabet
                std::back_inserter(sequence));
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

        detail::consume(stream_view | detail::take_line_or_throw);

        /* Qualities */
        auto qview = stream_view | std::views::filter(!is_space) // this consumes trailing newline
                   | detail::take_exactly_or_throw(sequence_size_after - sequence_size_before);
        if constexpr (!detail::decays_to_ignore_v<qual_type>)
        {
            std::ranges::copy(qview | views::char_to<std::ranges::range_value_t<qual_type>>,
                              std::back_inserter(qualities));
        }
        else
        {
            detail::consume(qview);
        }
#endif
    }

    //!\copydoc sequence_file_output_format::write_sequence_record
    template <typename stream_type, // constraints checked by file
              typename seq_type,    // other constraints checked inside function
              typename id_type,
              typename qual_type>
    void write_sequence_record(stream_type & stream,
                               sequence_file_output_options const & options,
                               seq_type && sequence,
                               id_type && id,
                               qual_type && qualities)
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
            throw std::logic_error{
                "The SEQ and SEQ_QUAL fields may not both be set to ignore when writing FASTQ files."};
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
            throw std::logic_error{
                "The QUAL and SEQ_QUAL fields may not both be set to ignore when writing FASTQ files."};
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

} // namespace seqan3
