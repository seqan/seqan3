// ============================================================================
//                 SeqAn - The Library for Sequence Analysis
// ============================================================================
//
// Copyright (c) 2006-2018, Knut Reinert & Freie Universitaet Berlin
// Copyright (c) 2016-2018, Knut Reinert & MPI Molekulare Genetik
// All rights reserved.
//
// Redistribution and use in source and binary forms, with or without
// modification, are permitted provided that the following conditions are met:
//
//     * Redistributions of source code must retain the above copyright
//       notice, this list of conditions and the following disclaimer.
//     * Redistributions in binary form must reproduce the above copyright
//       notice, this list of conditions and the following disclaimer in the
//       documentation and/or other materials provided with the distribution.
//     * Neither the name of Knut Reinert or the FU Berlin nor the names of
//       its contributors may be used to endorse or promote products derived
//       from this software without specific prior written permission.
//
// THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
// AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
// IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
// ARE DISCLAIMED. IN NO EVENT SHALL KNUT REINERT OR THE FU BERLIN BE LIABLE
// FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL
// DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR
// SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER
// CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT
// LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY
// OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH
// DAMAGE.
//
// ============================================================================

/*!\file
 * \brief Provides the seqan3::sequence_file_format_fasta class.
 * \author Hannes Hauswedell <hannes.hauswedell AT fu-berlin.de>
 */

#pragma once

#include <iterator>
#include <string>
#include <vector>

#include <range/v3/algorithm/copy.hpp>
#include <range/v3/utility/iterator.hpp>
#include <range/v3/view/drop_while.hpp>
#include <range/v3/view/remove_if.hpp>
#include <range/v3/view/transform.hpp>
#include <range/v3/view/take_while.hpp>

#include <seqan3/alphabet/nucleotide/dna5.hpp>
#include <seqan3/alphabet/quality/aliases.hpp>
#include <seqan3/core/metafunction/range.hpp>
#include <seqan3/io/detail/ignore_output_iterator.hpp>
#include <seqan3/io/detail/output_iterator_conversion_adaptor.hpp>
#include <seqan3/io/stream/parse_condition.hpp>
#include <seqan3/std/concept/range.hpp>
#include <seqan3/range/view/char_to.hpp>

namespace seqan3
{
/*!\brief       The FastA format.
 * \implements  sequence_file_format_concept
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
protected:
    //!\privatesection
    //!\brief Implementation of reading the ID.
    template <typename stream_view_t,
              typename options_type,
              typename id_type>
    void read_id(stream_view_t & stream_view,
                 options_type  & options,
                 id_type       & id)
    {
        auto const is_id = is_char<'>'>{} || is_char<';'>{};

        if (!is_id(*ranges::begin(stream_view)))
            throw parse_error{std::string{"Expected to be on beginning of ID, but "} + is_id.msg.string() +
                              " evaluated to false on " + detail::make_printable(*ranges::begin(stream_view))};

        // read id
        if (options.truncate_ids)
        {
            ranges::copy(stream_view | ranges::view::drop_while(is_id || is_blank)        // skip leading >
                                     | ranges::view::take_while(!(is_cntrl || is_blank)), // read ID until delimiter…
                         detail::make_conversion_output_iterator(id));                    // … ^A is old delimiter

            // consume rest of line
            ranges::copy(stream_view | ranges::view::take_while(!(is_char<'\n'>{} || is_char<'\r'>{})),
                         detail::make_conversion_output_iterator(std::ignore));
        }
        else
        {
            ranges::copy(stream_view | ranges::view::take_while(!(is_char<'\n'>{} || is_char<'\r'>{})) // read line
                                     | ranges::view::drop_while(is_id || is_blank),                    // skip leading >
                         detail::make_conversion_output_iterator(id));
        }
    }

    //!\brief Implementation of reading the sequence.
    template <typename      stream_view_t,
              typename      options_type,
              range_concept seq_type>
    void read_seq(stream_view_t  & stream_view,
                  options_type   &,
                  seq_type       & seq)
    {
        is_in_alphabet<typename options_type::sequence_legal_alphabet> const is_legal_alph;
        auto const is_id = is_char<'>'>{} || is_char<';'>{};

        ranges::copy(stream_view | ranges::view::take_while(!is_id)                // until next header (or end)
                                 | ranges::view::remove_if(is_space || is_digit)   // ignore whitespace and numbers
                                 | ranges::view::transform([is_legal_alph] (char const c)
                                   {
                                       if (!is_legal_alph(c))
                                       {
                                           throw parse_error{std::string{"Encountered an unexpected letter: "} +
                                                             is_legal_alph.msg.string() +
                                                             " evaluated to false on " +
                                                             detail::make_printable(c)};
                                       }
                                       return c;
                                   })                                              // enforce legal alphabet
                                 | view::char_to<value_type_t<seq_type>>,          // convert to actual target alphabet
                     detail::make_conversion_output_iterator(seq));
    }

    //!\brief Implementation of reading the sequence; overload for std::ignore (should be slightly faster).
    template <typename stream_view_t,
              typename options_type>
    void read_seq(stream_view_t          & stream_view,
                  options_type           &,
                  detail::ignore_t const &)
    {
        auto seq_view = stream_view | ranges::view::take_while(!is_char<'>'>{});     // until next header (or end)

        for (auto it = ranges::begin(seq_view); it != ranges::end(seq_view); ++it)
        {}
    }

public:
    /*!\name Constructors, destructor and assignment
     * \brief Rule of five explicitly defaulted.
     * \{
     */
    sequence_file_format_fasta() = default;
    sequence_file_format_fasta(sequence_file_format_fasta const &) = delete;
    sequence_file_format_fasta & operator=(sequence_file_format_fasta const &) = delete;
    sequence_file_format_fasta(sequence_file_format_fasta &&) = default;
    sequence_file_format_fasta & operator=(sequence_file_format_fasta &&) = default;
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

    //!\copydoc sequence_file_format_concept::read
    template <typename stream_type,     // constraints checked by file
              typename options_type,    // given by file
              typename seq_type,        // other constraints checked inside function
              typename id_type,
              typename qual_type,
              typename seq_qual_type>
    void read(stream_type        & stream,
              options_type const & options,
              seq_type           & sequence,
              id_type            & id,
              qual_type          & SEQAN3_DOXYGEN_ONLY(qualities),
              seq_qual_type      & seq_qual)
    {
        static_assert(detail::decays_to_ignore_v<seq_type> || detail::decays_to_ignore_v<seq_qual_type>,
                      "Either the sequence field, or the seq_qual field need to be set to std::ignore.");
        static_assert(detail::decays_to_ignore_v<qual_type> || detail::decays_to_ignore_v<seq_qual_type>,
                      "Either the qualities field, or the seq_qual field need to be set to std::ignore.");

        auto stream_view = ranges::iterator_range{std::istreambuf_iterator<char>{stream},
                                                  std::istreambuf_iterator<char>{}};
        // ID
        read_id(stream_view, options, id);

        // Sequence
        if constexpr (!detail::decays_to_ignore_v<seq_type>) // sequence
            read_seq(stream_view, options, sequence);
        else
            read_seq(stream_view, options, seq_qual);        // seq_qual (possibly std::ignore, too)

        // make sure "buffer at end" implies "stream at end"
        if ((std::istreambuf_iterator<char>{stream} == std::istreambuf_iterator<char>{}) &&
            (!stream.eof()))
        {
            stream.get(); // triggers error in stream and sets eof
        }
    }

    //!\copydoc sequence_file_format_concept::write
    template <typename stream_type,     // constraints checked by file
              typename options_type,    // given by file
              typename seq_type,        // other constraints checked inside function
              typename id_type,
              typename qual_type,
              typename seq_qual_type>
    void write([[maybe_unused]] stream_type        & stream,
               [[maybe_unused]] options_type const & options,
               [[maybe_unused]] seq_type           && sequence,
               [[maybe_unused]] id_type            && id,
               [[maybe_unused]] qual_type          && qualities,
               [[maybe_unused]] seq_qual_type      && seq_qual)
    {
        //TODO
    }
};

} // namespace seqan3

