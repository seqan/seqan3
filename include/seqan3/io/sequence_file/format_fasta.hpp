// -----------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2019, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2019, Knut Reinert & MPI für molekulare Genetik
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at: https://github.com/seqan/seqan3/blob/master/LICENSE.md
// -----------------------------------------------------------------------------------------------------

/*!\file
 * \brief Provides the seqan3::sequence_file_format_fasta class.
 * \author Hannes Hauswedell <hannes.hauswedell AT fu-berlin.de>
 */

#pragma once

#include <iterator>
#include <string>
#include <string_view>
#include <vector>

#include <range/v3/algorithm/copy.hpp>
#include <range/v3/view/chunk.hpp>
#include <range/v3/view/drop_while.hpp>
#include <range/v3/view/join.hpp>
#include <range/v3/view/remove_if.hpp>
#include <range/v3/view/take_while.hpp>

#include <seqan3/alphabet/nucleotide/dna5.hpp>
#include <seqan3/alphabet/quality/aliases.hpp>
#include <seqan3/core/metafunction/range.hpp>
#include <seqan3/io/detail/ignore_output_iterator.hpp>
#include <seqan3/io/detail/misc.hpp>
#include <seqan3/io/sequence_file/input_options.hpp>
#include <seqan3/io/sequence_file/output_options.hpp>
#include <seqan3/io/stream/parse_condition.hpp>
#include <seqan3/range/shortcuts.hpp>
#include <seqan3/range/detail/misc.hpp>
#include <seqan3/range/view/char_to.hpp>
#include <seqan3/range/view/to_char.hpp>
#include <seqan3/range/view/take.hpp>
#include <seqan3/range/view/take_exactly.hpp>
#include <seqan3/range/view/take_line.hpp>
#include <seqan3/range/view/take_until.hpp>
#include <seqan3/std/algorithm>
#include <seqan3/std/ranges>

namespace seqan3
{
/*!\brief       The FastA format.
 * \implements  SequenceFileInputFormat
 * \implements  SequenceFileOutputFormat
 * \ingroup     sequence
 *
 * \details
 *
 * ### Introduction
 *
 * FastA is the de-facto-standard for sequence storage in bionformatics. See the
 * [article on wikipedia](https://en.wikipedia.org/wiki/FASTA_format) for a an in-depth description of the format.
 *
 * ### Fields
 *
 * The FastA format provides the fields seqan3::field::SEQ and seqan3::field::ID. Both fields are required when writing.
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
class sequence_file_format_fasta
{
public:
    /*!\name Constructors, destructor and assignment
     * \{
     */
    sequence_file_format_fasta() = default;                                          //!< Defaulted
    //!\brief Copy construction is explicitly deleted, because you can't have multiple access to the same file.
    sequence_file_format_fasta(sequence_file_format_fasta const &) = delete;
    //!\brief Copy assignment is explicitly deleted, because you can't have multiple access to the same file.
    sequence_file_format_fasta & operator=(sequence_file_format_fasta const &) = delete;
    sequence_file_format_fasta(sequence_file_format_fasta &&) = default;             //!< Defaulted
    sequence_file_format_fasta & operator=(sequence_file_format_fasta &&) = default; //!< Defaulted
    ~sequence_file_format_fasta() = default;                                         //!< Defaulted
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

    //!\copydoc SequenceFileInputFormat::read
    template <typename stream_type,     // constraints checked by file
              typename seq_legal_alph_type, bool seq_qual_combined,
              typename seq_type,        // other constraints checked inside function
              typename id_type,
              typename qual_type>
    void read(stream_type                                                            & stream,
              sequence_file_input_options<seq_legal_alph_type, seq_qual_combined> const & options,
              seq_type                                                               & sequence,
              id_type                                                                & id,
              qual_type                                                              & SEQAN3_DOXYGEN_ONLY(qualities))
    {
        using stream_char_t = typename stream_type::char_type;
        auto stream_view = std::ranges::subrange<decltype(std::istreambuf_iterator<stream_char_t>{stream}),
                                          decltype(std::istreambuf_iterator<stream_char_t>{})>
                            {std::istreambuf_iterator<stream_char_t>{stream},
                             std::istreambuf_iterator<stream_char_t>{}};
        // ID
        read_id(stream_view, options, id);

        // Sequence
        read_seq(stream_view, options, sequence);
    }

    //!\copydoc SequenceFileOutputFormat::write
    template <typename stream_type,     // constraints checked by file
              typename seq_type,        // other constraints checked inside function
              typename id_type,
              typename qual_type>
    void write(stream_type                     & stream,
               sequence_file_output_options const & options,
               seq_type                       && sequence,
               id_type                        && id,
               qual_type                      && SEQAN3_DOXYGEN_ONLY(qualities))
    {

        std::ranges::ostreambuf_iterator stream_it{stream};

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

protected:
    //!\privatesection
    //!\brief Implementation of reading the ID.
    template <typename stream_view_t,
              typename seq_legal_alph_type, bool seq_qual_combined,
              typename id_type>
    void read_id(stream_view_t                                                          & stream_view,
                 sequence_file_input_options<seq_legal_alph_type, seq_qual_combined> const & options,
                 id_type                                                                & id)
    {
        auto const is_id = is_char<'>'> || is_char<';'>;

        if (!is_id(*begin(stream_view)))
            throw parse_error{std::string{"Expected to be on beginning of ID, but "} + is_id.msg.str() +
                              " evaluated to false on " + detail::make_printable(*begin(stream_view))};

        // read id
        if constexpr (!detail::decays_to_ignore_v<id_type>)
        {
            if (options.truncate_ids)
            {
                std::ranges::copy(stream_view | std::view::drop_while(is_id || is_blank)     // skip leading >
                                              | view::take_until_or_throw(is_cntrl || is_blank) // read ID until delimiter…
                                              | view::char_to<value_type_t<id_type>>,
                                  std::back_inserter(id));                                      // … ^A is old delimiter

                // consume rest of line
                detail::consume(stream_view | view::take_line_or_throw);
            }
            else
            {
                std::ranges::copy(stream_view | view::take_line_or_throw                                        // read line
                                              | std::view::drop_while(is_id || is_blank)                     // skip leading >
                                              | view::char_to<value_type_t<id_type>>,
                                  std::back_inserter(id));
            }
        }
        else
        {
            detail::consume(stream_view | view::take_line_or_throw);
        }
    }

    //!\brief Implementation of reading the sequence.
    template <typename stream_view_t,
              typename seq_legal_alph_type, bool seq_qual_combined,
              typename seq_type>
    void read_seq(stream_view_t                                                          & stream_view,
                  sequence_file_input_options<seq_legal_alph_type, seq_qual_combined> const &,
                  seq_type                                                               & seq)
    {
        auto constexpr is_id = is_char<'>'> || is_char<';'>;

        if constexpr (!detail::decays_to_ignore_v<seq_type>)
        {
            auto constexpr is_legal_alph = is_in_alphabet<seq_legal_alph_type>;
            std::ranges::copy(stream_view | view::take_until(is_id)                      // until next header (or end)
                                          | std::view::filter(!(is_space || is_digit))// ignore whitespace and numbers
                                          | std::view::transform([is_legal_alph] (char const c)
                                            {
                                                if (!is_legal_alph(c))
                                                {
                                                    throw parse_error{std::string{"Encountered an unexpected letter: "} +
                                                                        is_legal_alph.msg.str() +
                                                                        " evaluated to false on " +
                                                                        detail::make_printable(c)};
                                                }
                                                return c;
                                            })                                           // enforce legal alphabet
                                          | view::char_to<value_type_t<seq_type>>,       // convert to actual target alphabet
                              std::back_inserter(seq));
        }
        else
        {
            detail::consume(stream_view | view::take_until(is_id));
        }
    }

    //!\brief Implementation of writing the ID.
    template <typename stream_it_t,
              typename id_type>
    void write_id(stream_it_t                     & stream_it,
                  sequence_file_output_options const & options,
                  id_type                        && id)
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
    void write_seq(stream_it_t                    & stream_it,
                  sequence_file_output_options const & options,
                  seq_type                       && seq)
    {
        if (options.fasta_letters_per_line > 0)
        {
            std::ranges::copy(seq | view::to_char
                             | ranges::view::chunk(options.fasta_letters_per_line)
                             | std::view::join(options.add_carriage_return
                                                    ? std::string_view{"\r\n"}
                                                    : std::string_view{"\n"}),
                         stream_it);
            // TODO(h-2): benchmark the above vs:
//             size_t count = 0;
//             for (auto seq_it = begin(seq); seq_it != end(seq_it); ++seq_it)
//             {
//                 stream_it = to_char(*seq_it);
//                 ++count;
//                 if (count % fasta_letters_per_line == 0)
//                 {
//                     detail::write_eol(stream_it, options.add_carriage_return);
//                 }
//             }
        }
        else
        {
            std::ranges::copy(seq | view::to_char, stream_it);
        }

        detail::write_eol(stream_it, options.add_carriage_return);
    }
};

} // namespace seqan3
