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
 * \brief Provides the seqan3::sequence_file_format_sam class.
 * \author Mitra Darvish <mitra.darvish AT fu-berlin.de>
 */

#pragma once

#include <iterator>
#include <string>
#include <string_view>
#include <vector>

#include <range/v3/algorithm/copy.hpp>
#include <range/v3/utility/iterator.hpp>
#include <range/v3/view/chunk.hpp>
#include <range/v3/view/join.hpp>
#include <range/v3/view/remove_if.hpp>

#include <seqan3/alphabet/nucleotide/dna5.hpp>
#include <seqan3/alphabet/quality/aliases.hpp>
#include <seqan3/core/metafunction/range.hpp>
#include <seqan3/io/detail/ignore_output_iterator.hpp>
#include <seqan3/io/detail/misc.hpp>
#include <seqan3/io/sequence_file/input_options.hpp>
#include <seqan3/io/sequence_file/output_options.hpp>
#include <seqan3/io/stream/parse_condition.hpp>
#include <seqan3/range/detail/misc.hpp>
#include <seqan3/range/view/char_to.hpp>
#include <seqan3/range/view/to_char.hpp>
#include <seqan3/range/view/take.hpp>
#include <seqan3/range/view/take_exactly.hpp>
#include <seqan3/range/view/take_line.hpp>
#include <seqan3/range/view/take_until.hpp>
#include <seqan3/std/ranges>
#include <seqan3/std/view/subrange.hpp>
#include <seqan3/std/view/transform.hpp>

namespace seqan3
{
  /*!\brief       The sam format used as sequence file.
   * \implements  sequence_file_format_concept
   * \ingroup     sequence
   *
   * \details
   *
   * ### Introduction
   *
   * The SAM format is commonly used to store pairwise alignment information between a query sequence and its
   * reference sequence, e.g. a read mapping result. As it provides information on the **id**, **sequence** and
   * **quality** of each query, it can additionally be treated as a sequence file.
   * See the [article on wikipedia](https://en.wikipedia.org/wiki/SAM_(file_format)) for an in-depth description of
   * the format.
   *
   * ### Fields
   *
   * The SAM format provides the fields seqan3::field::SEQ, seqan3::field::ID and seqan3::field::SEQ_QUAL.
   * All fields are allowed to be empty when writing.
   *
   * ### Implementation notes
   *
   * This implementation supports the following optional features of the format:
   *
   * line breaks and/or other whitespace characters in any part of the sequence (only when reading!), optional tags
   * are going to be ignored (when reading)
   *
   */
class sequence_file_format_sam
{
public:
    /*!\name Constructors, destructor and assignment
     * \{
     */
     //!\brief Default constructor is explicitly deleted, you need to give a stream or file name.
    sequence_file_format_sam() = default;
    //!\brief Copy construction is explicitly deleted, because you can't have multiple access to the same file.
    sequence_file_format_sam(sequence_file_format_sam const &) = delete;
    //!\brief Copy assignment is explicitly deleted, because you can't have multiple access to the same file.
    sequence_file_format_sam & operator=(sequence_file_format_sam const &) = delete;
    //!\brief Move construction is defaulted.
    sequence_file_format_sam(sequence_file_format_sam &&) = default;
    //!\brief Move assignment is defaulted.
    sequence_file_format_sam & operator=(sequence_file_format_sam &&) = default;
    //!\}

    //!\brief The valid file extensions for this format; note that you can modify this value.
    static inline std::vector<std::string> file_extensions
    {
        { "sam" },
    };

    //!\copydoc sequence_file_input_format_concept::read
    template <typename stream_type,     // constraints checked by file
              typename seq_legal_alph_type, bool seq_qual_combined,
              typename seq_type,        // other constraints checked inside function
              typename id_type,
              typename qual_type>
    void read(stream_type                     & stream,
              sequence_file_input_options<seq_legal_alph_type, seq_qual_combined> const & options,
              seq_type                        & sequence,
              id_type                         & id,
              qual_type                       & qualities)
    {
        auto stream_view = view::subrange<decltype(std::istreambuf_iterator<char>{stream}),
                                          decltype(std::istreambuf_iterator<char>{})>
                            {std::istreambuf_iterator<char>{stream},
                             std::istreambuf_iterator<char>{}};

        // cache the begin position so we write quals to the same position as seq in seq_qual case
        size_t sequence_size_before = 0;
        size_t sequence_size_after = 0;
        if constexpr (!detail::decays_to_ignore_v<seq_type>)
            sequence_size_before = ranges::size(sequence);

        auto constexpr is_tab = is_char<'\t'>;
        auto constexpr is_cntrl= is_char<'\n'>;

        while (is_char<'@'>(*ranges::begin(stream_view)))
            detail::consume(stream_view | view::take_line);

        /* ID */
        if constexpr (!detail::decays_to_ignore_v<id_type>)
        {
            if(is_char<'*'>(*begin(stream_view)))
                throw parse_error{"The ID field may not be empty for SAM files."};

            if (options.truncate_ids)
            {
              ranges::copy(stream_view | view::take_until_or_throw(is_blank)
                                       | view::char_to<value_type_t<id_type>>,
                           std::back_inserter(id));
              detail::consume(stream_view | view::take_until_or_throw(is_tab));
            }
            else
            {
              ranges::copy(stream_view | view::take_until_or_throw(is_tab)
                                       | view::char_to<value_type_t<id_type>>,
                           std::back_inserter(id));
            }
        }
        else
        {
             detail::consume(stream_view | view::take_until_or_throw(is_tab));
        }

        //Jump over SAM colums 2-9 (FLAG RNAME POS MAPQ CIGAR RNEXT PNEXT TLEN)
        for (int i = 0; i < 8;i++)
        {
            ++begin(stream_view);
            detail::consume(stream_view | view::take_until_or_throw(is_tab));
        }
        ++begin(stream_view);

        auto seq_view = stream_view | view::take_until_or_throw(is_tab); // until next tab

         // Sequence
         if constexpr (!detail::decays_to_ignore_v<seq_type>)
         {
             if(is_char<'*'>(*begin(stream_view)))
                 throw parse_error{"The Sequence field may not be empty for SAM files."};

             auto constexpr is_legal_alph = is_in_alphabet<seq_legal_alph_type>;
             ranges::copy(seq_view | view::transform([is_legal_alph] (char const c) // enforce legal alphabet
                                      {
                                          if (!is_legal_alph(c))
                                          {
                                              throw parse_error{std::string{"Encountered an unexpected letter: "} +
                                                                is_legal_alph.msg.string() +
                                                                " evaluated to false on " +
                                                                detail::make_printable(c)};
                                          }
                                          return c;
                                      })
                                    | view::char_to<value_type_t<seq_type>>, // convert to actual target alphabet
                           std::back_inserter(sequence));

              sequence_size_after = ranges::size(sequence);
          }
          else
          {
              for (auto it = begin(seq_view); it != end(seq_view); ++it)
                  ++sequence_size_after;
          }
          ++begin(stream_view); // skip tab

          /* Qualities */
          if (is_char<'*'>(*begin(stream_view)))
          {
              ++begin(stream_view); // skip *
          }
          else
          {
             auto qual = stream_view | view::take_until_or_throw(is_tab || is_cntrl);
             auto qual_view = qual | view::take_exactly_or_throw(sequence_size_after - sequence_size_before);

             if constexpr (seq_qual_combined)
             {
                  // seq_qual field implies that they are the same variable
                  assert(std::addressof(sequence) == std::addressof(qualities));
                  ranges::copy(qual_view | view::char_to<typename value_type_t<qual_type>::quality_alphabet_type>,
                               ranges::begin(qualities) + sequence_size_before);
             }
             else if constexpr (!detail::decays_to_ignore_v<qual_type>)
             {
                  ranges::copy(qual_view | view::char_to<value_type_t<qual_type>>,
                             std::back_inserter(qualities));
             }
             else
             {
                  detail::consume(qual_view);
             }

             if (!(is_tab(*begin(qual)) | is_cntrl(*begin(qual))))
                throw unexpected_end_of_input{"Quality length surpasses sequence length."};
          }

          // consume the remaining characters (optional tags and newline)
          detail::consume(stream_view | view::take_until_or_throw(is_cntrl));
          ++begin(stream_view);

          // make sure "buffer at end" implies "stream at end"
          if ((std::istreambuf_iterator<char>{stream} == std::istreambuf_iterator<char>{}) &&
              (!stream.eof()))
          {
              stream.get(); // triggers error in stream and sets eof
          }
    }

    //!\copydoc sequence_file_output_format_concept::write
    template <typename stream_type,     // constraints checked by file
              typename seq_type,        // other constraints checked inside function
              typename id_type,
              typename qual_type>
    void write(stream_type                     & stream,
               sequence_file_output_options const & options,
               seq_type                       && sequence,
               id_type                        && id,
               qual_type                      && qualities)
    {
        if constexpr (!(detail::decays_to_ignore_v<seq_type>))
        {
            static_assert(std::ranges::ForwardRange<seq_type> && alphabet_concept<value_type_t<seq_type>>,
                          "The sequence must model std::ranges::ForwardRange and its value type must model"
                          " seqan3::alphabet_concept.");
        }
        ranges::ostreambuf_iterator stream_it{stream};
        // ID
        if constexpr (detail::decays_to_ignore_v<id_type>)
        {
            stream_it = '*';
        }
        else
        {
            if (ranges::empty(id)) //[[unlikely]]
                stream_it = '*';
            ranges::copy(id, stream_it);
        }

        stream << "\t0\t*\t0\t0\t*\t*\t0\t0\t";

        // Sequence
        if constexpr (detail::decays_to_ignore_v<seq_type>)
        {
            stream_it = '*';
        }
        else
        {
            if (ranges::empty(sequence)) //[[unlikely]]
                stream_it = '*';
            ranges::copy(sequence | view::to_char, stream_it);
        }
        stream_it = '\t';

        // Quality line
        if constexpr (detail::decays_to_ignore_v<qual_type>)
          stream_it = '*';
        else if (ranges::empty(qualities))
          stream_it = '*';
        else
            ranges::copy(qualities | view::to_char, stream_it);
        detail::write_eol(stream_it, options.add_carriage_return);
    }
};

} // namespace seqan3
