// -----------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2019, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2019, Knut Reinert & MPI für molekulare Genetik
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at: https://github.com/seqan/seqan3/blob/master/LICENSE.md
// -----------------------------------------------------------------------------------------------------

/*!\file
 * \brief Provides the seqan3::format_fastq tag and the seqan3::sequence_file_input_format and
 *        seqan3::sequence_file_output_format specialisation for this tag.
 * \author Hannes Hauswedell <hannes.hauswedell AT fu-berlin.de>
 */

#pragma once

#include <iterator>
#include <string>
#include <string_view>
#include <vector>

#include <range/v3/algorithm/copy.hpp>
#include <range/v3/view/chunk.hpp>
#include <range/v3/view/join.hpp>
#include <range/v3/view/remove_if.hpp>

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
#include <seqan3/range/shortcuts.hpp>
#include <seqan3/range/detail/misc.hpp>
#include <seqan3/range/view/char_to.hpp>
#include <seqan3/range/view/istreambuf.hpp>
#include <seqan3/range/view/to_char.hpp>
#include <seqan3/range/view/take.hpp>
#include <seqan3/range/view/take_exactly.hpp>
#include <seqan3/range/view/take_line.hpp>
#include <seqan3/range/view/take_until.hpp>
#include <seqan3/std/algorithm>
#include <seqan3/std/ranges>

namespace seqan3
{

/*!\brief       The FastQ format. (tag)
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
 * ### Fields
 *
 * The FastQ format provides the fields seqan3::field::SEQ, seqan3::field::ID and seqan3::field::QUAL; or alternatively
 * provides seqan3::field::SEQ_QUAL as a single field of sequence and quality. All three fields (or ID + SEQ_QUAL) are
 * required when writing and the sequence and qualities are required to be of the same length.
 *
 * ### Encodings
 *
 * All documented encodings for the quality string are supported (see the article above), but they are **not detected**
 * from the file. Instead, when reading the file, you have to set the respective alphabet via a traits type (see
 * seqan3::SequenceFileInputTraits and the \ref quality submodule).
 *
 * ### Implementation notes
 *
 * This implementation supports the following optional features of the format:
 *
 *   * line breaks and/or other whitespace characters in any part of the sequence and/or qualities (only when reading!)
 *   * writing the ID to the `+`-line also (line is always ignored when reading)
 *
 */
struct format_fastq
{
    //!\brief The valid file extensions for this format; note that you can modify this value.
    static inline std::vector<std::string> file_extensions
    {
        { "fastq" },
        { "fq"    }
    };
};

} // namespace seqan

namespace seqan3::detail
{

/*!\brief The seqan3::sequence_file_input_format specialisation that handles formatted FASTQ input.
 * \ingroup sequence
 * \tparam stream_char_type The underlying character type of the stream (usually `char`).
 */
template <typename stream_char_type>
class sequence_file_input_format<format_fastq, stream_char_type>
{
public:
    //!\brief Exposes the format tag that this class is specialised with.
    using format_tag = format_fastq;

    /*!\name Constructors, destructor and assignment
     * \{
     */
    sequence_file_input_format()                                                        = default; //!< Defaulted.
    //!\brief Copy construction is explicitly deleted, because you can't have multiple access to the same file.
    sequence_file_input_format(sequence_file_input_format const &)                      = delete;
    //!\brief Copy assignment is explicitly deleted, because you can't have multiple access to the same file.
    sequence_file_input_format & operator=(sequence_file_input_format const &)          = delete;
    sequence_file_input_format(sequence_file_input_format &&)                  noexcept = default; //!< Defaulted.
    sequence_file_input_format & operator=(sequence_file_input_format &&)      noexcept = default; //!< Defaulted.
    ~sequence_file_input_format()                                              noexcept = default; //!< Defaulted.

    //!\brief Construct the format given an input stream to read from.
    sequence_file_input_format(std::basic_istream<stream_char_type> & stream) :
        stream_view{view::istreambuf(stream)}
    {}
    //!\}

    //!\copydoc SequenceFileInputFormat::read
    template <typename seq_legal_alph_type, bool seq_qual_combined,
              typename seq_type,        // constraints checked by file, other constraints checked inside function
              typename id_type,
              typename qual_type>
    bool read(sequence_file_input_options<seq_legal_alph_type, seq_qual_combined> const & options,
              seq_type                                                                  & sequence,
              id_type                                                                   & id,
              qual_type                                                                 & qualities)
    {
        if (std::ranges::begin(stream_view) == std::ranges::end(stream_view)) // file has no records any more
            return true;

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
                std::ranges::copy(stream_view | view::take_until_or_throw(is_cntrl || is_blank)
                                              | view::char_to<value_type_t<id_type>>,
                                  std::back_inserter(id));
                detail::consume(stream_view | view::take_line_or_throw);
            }
            else
            {
                std::ranges::copy(stream_view | view::take_line_or_throw
                                              | view::char_to<value_type_t<id_type>>,
                                  std::back_inserter(id));
            }
        }
        else
        {
            detail::consume(stream_view | view::take_line_or_throw);
        }

        /* Sequence */
        auto seq_view = stream_view | view::take_until_or_throw(is_char<'+'>)    // until 2nd ID line
                                    | std::view::filter(!is_space);           // ignore whitespace
        if constexpr (!detail::decays_to_ignore_v<seq_type>)
        {
            auto constexpr is_legal_alph = is_in_alphabet<seq_legal_alph_type>;
            std::ranges::copy(seq_view | std::view::transform([is_legal_alph] (char const c) // enforce legal alphabet
                                    {
                                        if (!is_legal_alph(c))
                                        {
                                            throw parse_error{std::string{"Encountered an unexpected letter: "} +
                                                                is_legal_alph.msg.str() +
                                                                " evaluated to false on " +
                                                                detail::make_printable(c)};
                                        }
                                        return c;
                                    })
                                        | view::char_to<value_type_t<seq_type>>,         // convert to actual target alphabet
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

        /* 2nd ID line */
        if (*stream_it != '+') // [[unlikely]]
        {
            throw parse_error{std::string{"Expected '+' on beginning of 2nd ID line, got: "} +
                              detail::make_printable(*stream_it)};
        }
        detail::consume(stream_view | view::take_line_or_throw);

        /* Qualities */
        auto qview = stream_view | std::view::filter(!is_space)                  // this consumes trailing newline
                                 | view::take_exactly_or_throw(sequence_size_after - sequence_size_before);
        if constexpr (seq_qual_combined)
        {
            // seq_qual field implies that they are the same variable
            assert(std::addressof(sequence) == std::addressof(qualities));
            std::ranges::copy(qview | view::char_to<typename value_type_t<qual_type>::quality_alphabet_type>,
                              begin(qualities) + sequence_size_before);
        }
        else if constexpr (!detail::decays_to_ignore_v<qual_type>)
        {
            std::ranges::copy(qview | view::char_to<value_type_t<qual_type>>,
                              std::back_inserter(qualities));
        }
        else
        {
            detail::consume(qview);
        }

        return false;
    }
private:
    //!\brief A view over the file stream.
    decltype(view::istreambuf(std::declval<std::basic_istream<stream_char_type> &>())) stream_view{};
};

/*!\brief The seqan3::sequence_file_output_format specialisation that can write formatted FASTQ.
 * \ingroup sequence
 * \tparam stream_char_type The underlying character type of the stream (usually `char`).
 */
template <typename stream_char_type>
class sequence_file_output_format<format_fastq, stream_char_type>
{
public:
    //!\brief Exposes the format tag that this class is specialised with.
    using format_tag = format_fastq;

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

    //!\brief Construct from an output stream to write to.
    sequence_file_output_format(std::basic_ostream<stream_char_type> & stream) :
        stream_it{stream}
    {}
    //!\}

    //!\copydoc SequenceFileOutputFormat::write
    template <typename seq_type,        // constraints checked by file, other constraints checked inside function
              typename id_type,
              typename qual_type>
    void write(sequence_file_output_options const & options,
               seq_type                          && sequence,
               id_type                           && id,
               qual_type                         && qualities)
    {
        // ID
        if constexpr (detail::decays_to_ignore_v<id_type>)
        {
            throw std::logic_error{"The ID field may not be set to ignore when writing FASTQ files."};
        }
        else
        {
            if (empty(id)) //[[unlikely]]
                throw std::runtime_error{"The ID field may not be empty when writing FASTQ files."};

            stream_it = '@';
            std::ranges::copy(id, stream_it);

            detail::write_eol(stream_it, options.add_carriage_return);
        }

        // Sequence
        if constexpr (detail::decays_to_ignore_v<seq_type>)
        {
            throw std::logic_error{"The SEQ and SEQ_QUAL fields may not both be set to ignore when writing FASTQ files."};
        }
        else
        {
            if (empty(sequence)) //[[unlikely]]
                throw std::runtime_error{"The SEQ field may not be empty when writing FASTQ files."};

            std::ranges::copy(sequence | view::to_char, stream_it);

            detail::write_eol(stream_it, options.add_carriage_return);
        }

        // 2nd ID-line
        if constexpr (!detail::decays_to_ignore_v<id_type>)
        {
            stream_it = '+';

            if (options.fastq_double_id)
                std::ranges::copy(id, stream_it);

            detail::write_eol(stream_it, options.add_carriage_return);
        }

        // Quality line
        if constexpr (detail::decays_to_ignore_v<qual_type>)
        {
            throw std::logic_error{"The QUAL and SEQ_QUAL fields may not both be set to ignore when writing FASTQ files."};
        }
        else
        {
            if (empty(qualities)) //[[unlikely]]
                throw std::runtime_error{"The SEQ field may not be empty when writing FASTQ files."};

            if constexpr (std::ranges::SizedRange<seq_type> && std::ranges::SizedRange<qual_type>)
            {
                assert(size(sequence) == size(qualities));
            }

            std::ranges::copy(qualities | view::to_char, stream_it);

            detail::write_eol(stream_it, options.add_carriage_return);
        }
    }

private:
    //!\brief An ostreambuf iterator to the output stream to write to.
    seqan3::ostreambuf_iterator<stream_char_type> stream_it{};
};

} // namespace seqan3::detail
