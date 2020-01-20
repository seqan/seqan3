// -----------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2020, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2020, Knut Reinert & MPI für molekulare Genetik
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at: https://github.com/seqan/seqan3/blob/master/LICENSE.md
// -----------------------------------------------------------------------------------------------------

/*!\file
 * brief Provides the seqan3::format_fasta.
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
#include <seqan3/range/shortcuts.hpp>
#include <seqan3/range/detail/misc.hpp>
#include <seqan3/range/views/char_to.hpp>
#include <seqan3/range/views/istreambuf.hpp>
#include <seqan3/range/views/join.hpp>
#include <seqan3/range/views/to_char.hpp>
#include <seqan3/range/views/take.hpp>
#include <seqan3/range/views/take_exactly.hpp>
#include <seqan3/range/views/take_line.hpp>
#include <seqan3/range/views/take_until.hpp>
#include <seqan3/std/algorithm>
#include <seqan3/std/ranges>

namespace seqan3
{

/*!\brief       The FastA format.
 * \implements  sequence_file_input_format
 * \implements  sequence_file_output_format
 * \ingroup     sequence
 *
 * \details
 *
 * ### Introduction
 *
 * FastA is the de-facto-standard for sequence storage in bionformatics. See the
 * [article on wikipedia](https://en.wikipedia.org/wiki/FASTA_format) for a an in-depth description of the format.
 *
 * ### fields_specialisation
 *
 * The FastA format provides the fields seqan3::field::seq and seqan3::field::id. Both fields are required when writing.
 *
 * ### Implementation notes
 *
 * When reading the ID-line the identifier (either `;` or `>`) and any blank characters before the actual ID are
 * stripped.
 *
 * This implementation supports the following less known and optional features of the format:
 *
 *   * ID lines beginning with `;` instead of `>`
 *   * line breaks and other whitespace characters in any part of the sequence
 *   * character counts within the sequence (they are simply ignored)
 *
 * The following optional features are currently **not supported:**
 *
 *   * Multiple comment lines (starting with either `;` or `>`), only one ID line before the sequence line is accepted
 *
 */
class format_fasta
{
public:
    /*!\name Constructors, destructor and assignment
     * \{
     */
    format_fasta() noexcept = default; //!< Defaulted.
    format_fasta(format_fasta const &) noexcept = default; //!< Defaulted.
    format_fasta & operator=(format_fasta const &) noexcept = default; //!< Defaulted.
    format_fasta(format_fasta &&) noexcept = default; //!< Defaulted.
    format_fasta & operator=(format_fasta &&) noexcept = default; //!< Defaulted.
    ~format_fasta() noexcept = default; //!< Defaulted.
    //!\}

    //!\brief The valid file extensions for this format; note that you can modify this value.
    static inline std::vector<std::string> file_extensions
    {
        { "fasta" },
        { "fa"    },
        { "fna"   },
        { "ffn"   },
        { "faa"   },
        { "frn"   },
    };

protected:
    //!\copydoc sequence_file_input_format::read_sequence_record
    template <typename stream_type,     // constraints checked by file
              typename legal_alph_type, bool seq_qual_combined,
              typename seq_type,        // other constraints checked inside function
              typename id_type,
              typename qual_type>
    void read_sequence_record(stream_type & stream,
                              sequence_file_input_options<legal_alph_type, seq_qual_combined> const & options,
                              seq_type & sequence,
                              id_type & id,
                              qual_type & SEQAN3_DOXYGEN_ONLY(qualities))
    {
        auto stream_view = views::istreambuf(stream);

        // ID
        read_id(stream_view, options, id);

        // Sequence
        read_seq(stream_view, options, sequence);
    }

    //!\copydoc sequence_file_output_format::write_sequence_record
    template <typename stream_type,     // constraints checked by file
              typename seq_type,        // other constraints checked inside function
              typename id_type,
              typename qual_type>
    void write_sequence_record(stream_type & stream,
                               sequence_file_output_options const & options,
                               seq_type && sequence,
                               id_type && id,
                               qual_type && SEQAN3_DOXYGEN_ONLY(qualities))
    {
        seqan3::ostreambuf_iterator stream_it{stream};

        // ID
        if constexpr (detail::decays_to_ignore_v<id_type>)
        {
            throw std::logic_error{"The ID field may not be set to ignore when writing FASTA files."};
        }
        else
        {
            if (empty(id)) //[[unlikely]]
                throw std::runtime_error{"The ID field may not be empty when writing FASTA files."};

            write_id(stream_it, options, id);
        }

        // Sequence
        if constexpr (detail::decays_to_ignore_v<seq_type>) // sequence
        {
            throw std::logic_error{"The SEQ and SEQ_QUAL fields may not both be set to ignore when writing FASTA files."};
        }
        else
        {
            if (empty(sequence)) //[[unlikely]]
                throw std::runtime_error{"The SEQ field may not be empty when writing FASTA files."};

            write_seq(stream_it, options, sequence);
        }
    }

private:
    //!\privatesection
    //!\brief Implementation of reading the ID.
    template <typename stream_view_t,
              typename seq_legal_alph_type, bool seq_qual_combined,
              typename id_type>
    void read_id(stream_view_t & stream_view,
                  sequence_file_input_options<seq_legal_alph_type, seq_qual_combined> const & options,
                  id_type & id)
    {
        auto const is_id = is_char<'>'> || is_char<';'>;

        if (!is_id(*begin(stream_view)))
            throw parse_error{std::string{"Expected to be on beginning of ID, but "} + is_id.msg +
                              " evaluated to false on " + detail::make_printable(*begin(stream_view))};

        // read id
        if constexpr (!detail::decays_to_ignore_v<id_type>)
        {
            if (options.truncate_ids)
            {
            #if SEQAN3_WORKAROUND_VIEW_PERFORMANCE
                auto it = stream_view.begin();
                auto e = stream_view.end();
                for (; (it != e) && (is_id || is_blank)(*it); ++it)
                {}

                bool at_delimiter = false;
                for (; it != e; ++it)
                {
                    if ((is_cntrl || is_blank)(*it))
                    {
                        at_delimiter = true;
                        break;
                    }
                    id.push_back(assign_char_to(*it, value_type_t<id_type>{}));
                }

                if (!at_delimiter)
                    throw unexpected_end_of_input{"FastA ID line did not end in newline."};

                for (; (it != e) && ((!is_char<'\n'>)(*it)); ++it)
                {}

            #else // ↑↑↑ WORKAROUND | ORIGINAL ↓↓↓

                std::ranges::copy(stream_view | std::views::drop_while(is_id || is_blank)        // skip leading >
                                              | views::take_until_or_throw(is_cntrl || is_blank) // read ID until delimiter…
                                              | views::char_to<value_type_t<id_type>>,
                                  std::ranges::back_inserter(id));                               // … ^A is old delimiter

                // consume rest of line
                detail::consume(stream_view | views::take_line_or_throw);
            #endif // SEQAN3_WORKAROUND_VIEW_PERFORMANCE

            }
            else
            {
            #if SEQAN3_WORKAROUND_VIEW_PERFORMANCE
                auto it = stream_view.begin();
                auto e = stream_view.end();
                for (; (it != e) && (is_id || is_blank)(*it); ++it)
                {}

                bool at_delimiter = false;
                for (; it != e; ++it)
                {
                    if ((is_char<'\n'>)(*it))
                    {
                        at_delimiter = true;
                        break;
                    }
                    id.push_back(assign_char_to(*it, value_type_t<id_type>{}));
                }

                if (!at_delimiter)
                    throw unexpected_end_of_input{"FastA ID line did not end in newline."};

            #else // ↑↑↑ WORKAROUND | ORIGINAL ↓↓↓

                std::ranges::copy(stream_view | views::take_line_or_throw                    // read line
                                              | std::views::drop_while(is_id || is_blank)    // skip leading >
                                              | views::char_to<value_type_t<id_type>>,
                                  std::ranges::back_inserter(id));
            #endif // SEQAN3_WORKAROUND_VIEW_PERFORMANCE
            }
        }
        else
        {
            detail::consume(stream_view | views::take_line_or_throw);
        }
    }

    //!\brief Implementation of reading the sequence.
    template <typename stream_view_t,
              typename seq_legal_alph_type, bool seq_qual_combined,
              typename seq_type>
    void read_seq(stream_view_t & stream_view,
                   sequence_file_input_options<seq_legal_alph_type, seq_qual_combined> const &,
                   seq_type & seq)
    {
        auto constexpr is_id = is_char<'>'> || is_char<';'>;

        if constexpr (!detail::decays_to_ignore_v<seq_type>)
        {
            auto constexpr not_in_alph = !is_in_alphabet<seq_legal_alph_type>;

        #if SEQAN3_WORKAROUND_VIEW_PERFORMANCE
            auto it = stream_view.begin();
            auto e = stream_view.end();
            for (; (it != e) && ((!is_id)(*it)); ++it)
            {
                if ((is_space || is_digit)(*it))
                    continue;
                else if (not_in_alph(*it))
                {
                    throw parse_error{std::string{"Encountered an unexpected letter: "} +
                                        not_in_alph.msg +
                                        " evaluated to true on " +
                                        detail::make_printable(*it)};
                }

                seq.push_back(assign_char_to(*it, value_type_t<seq_type>{}));
            }

        #else // ↑↑↑ WORKAROUND | ORIGINAL ↓↓↓

            std::ranges::copy(stream_view | views::take_until(is_id)                   // until next header (or end)
                                          | std::views::filter(!(is_space || is_digit))// ignore whitespace and numbers
                                          | std::views::transform([not_in_alph] (char const c)
                                            {
                                                if (not_in_alph(c))
                                                {
                                                    throw parse_error{std::string{"Encountered an unexpected letter: "} +
                                                                        not_in_alph.msg +
                                                                        " evaluated to false on " +
                                                                        detail::make_printable(c)};
                                                }
                                                return c;
                                            })                                      // enforce legal alphabet
                                          | views::char_to<value_type_t<seq_type>>, // convert to actual target alphabet
                              std::ranges::back_inserter(seq));
        #endif // SEQAN3_WORKAROUND_VIEW_PERFORMANCE
        }
        else
        {
            detail::consume(stream_view | views::take_until(is_id));
        }
    }

    //!\brief Implementation of writing the ID.
    template <typename stream_it_t,
              typename id_type>
    void write_id(stream_it_t & stream_it,
                   sequence_file_output_options const & options,
                   id_type && id)
    {
        if (options.fasta_legacy_id_marker)
            stream_it = ';';
        else
            stream_it = '>';

        if (options.fasta_blank_before_id)
            stream_it = ' ';

        std::ranges::copy(id, stream_it);

        detail::write_eol(stream_it, options.add_carriage_return);
    }

    //!\brief Implementation of writing the sequence.
    template <typename stream_it_t,
              typename seq_type>
    void write_seq(stream_it_t & stream_it,
                   sequence_file_output_options const & options,
                   seq_type && seq)
    {
        if (options.fasta_letters_per_line > 0)
        {
        #if SEQAN3_WORKAROUND_VIEW_PERFORMANCE
            size_t count = 0;
            for (auto c : seq)
            {
                stream_it = to_char(c);
                if (++count % options.fasta_letters_per_line == 0)
                    detail::write_eol(stream_it, options.add_carriage_return);
            }
            if (count % options.fasta_letters_per_line != 0)
                detail::write_eol(stream_it, options.add_carriage_return);

        #else // ↑↑↑ WORKAROUND | ORIGINAL ↓↓↓

            //TODO: combining chunk and join is substantially faster than views::interleave (2.5x), why?
            std::ranges::copy(seq | views::to_char
                                  | ranges::view::chunk(options.fasta_letters_per_line)
                                  | views::join(options.add_carriage_return
                                                    ? std::string_view{"\r\n"}
                                                    : std::string_view{"\n"}),
                              stream_it);
            detail::write_eol(stream_it, options.add_carriage_return);
        #endif // SEQAN3_WORKAROUND_VIEW_PERFORMANCE
        }
        else
        {
            // No performance workaround here, because transform views alone are fast
            std::ranges::copy(seq | views::to_char, stream_it);
            detail::write_eol(stream_it, options.add_carriage_return);
        }
    }
};

} // namespace seqan
