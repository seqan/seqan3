// -----------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2019, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2019, Knut Reinert & MPI für molekulare Genetik
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at: https://github.com/seqan/seqan3/blob/master/LICENSE
// -----------------------------------------------------------------------------------------------------

/*!\file
 * \brief Provides the seqan3::sequence_file_format_sam class.
 * \author Mitra Darvish <mitra.darvish AT fu-berlin.de>
 */

#pragma once

#include <iterator>
#include <string>
#include <string_view>
#include <vector>

#include <seqan3/alphabet/nucleotide/dna5.hpp>
#include <seqan3/core/metafunction/range.hpp>
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

namespace seqan3
{
  /*!\brief       The SAM format used as sequence file.
   * \implements  SequenceFileFormat
   * \ingroup     sequence
   *
   * \details
   *
   * ### Introduction
   *
   * The SAM format is commonly used to store pairwise alignment information between a query sequence and its
   * reference sequence, e.g. a read mapping result. Some people also use the SAM format as plain storage for
   * sequences (and qualities) and in some cases the original sequence files are no longer available. The
   * seqan3::sequence_file_format_sam allows using SAM files in this manner and provides easy convertibility
   * from/to FASTQ; but there is no access to the alignment information stored in SAM files. Use
   * seqan3::alignment_file_format_sam if you are interested in the alignment.
   * See the [article on Wikipedia](https://en.wikipedia.org/wiki/SAM_(file_format)) or the
   * [technical specification] (https://samtools.github.io/hts-specs/SAMv1.pdf) for an in-depth description of
   * the format.
   *
   * ### Fields
   *
   * The SAM format provides the fields seqan3::field::SEQ, seqan3::field::ID and seqan3::field::SEQ_QUAL.
   * All fields are allowed to be empty when writing.
   *
   * ### Implementation notes
   *
   * This implementation ignores all fields besides id, seq and quality.
   *
   */
class sequence_file_format_sam
{
public:
    /*!\name Constructors, destructor and assignment
     * \{
     */
     //!\brief Default constructor is explicitly defaulted, you need to give a stream or file name.
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

    //!\copydoc SequenceFileInputFormat::read
    template <typename stream_type,     // constraints checked by file
              typename seq_legal_alph_type, bool seq_qual_combined,
              typename seq_type,        // other constraints checked inside function
              typename id_type,
              typename qual_type>
    void read(stream_type                                                               & stream,
              sequence_file_input_options<seq_legal_alph_type, seq_qual_combined> const & options,
              seq_type                                                                  & sequence,
              id_type                                                                   & id,
              qual_type                                                                 & qualities)
    {
        auto stream_view = std::ranges::subrange<decltype(std::istreambuf_iterator<char>{stream}),
                                          decltype(std::istreambuf_iterator<char>{})>
                            {std::istreambuf_iterator<char>{stream},
                             std::istreambuf_iterator<char>{}};

        // cache the begin position so we write quals to the same position as seq in seq_qual case
        size_t sequence_size_before{0};
        size_t sequence_size_after{0};
        if constexpr (!detail::decays_to_ignore_v<seq_type>)
            sequence_size_before = std::ranges::size(sequence);

        auto constexpr is_tab{is_char<'\t'>};
        auto constexpr is_cntrl{is_char<'\n'>};

        while (is_char<'@'>(*ranges::begin(stream_view)))
            detail::consume(stream_view | view::take_line);

        /* ID */
        if (is_char<'*'>(*begin(stream_view)))
            throw parse_error{"The ID field may not be empty for sequence files."};

        if constexpr (!detail::decays_to_ignore_v<id_type>)
        {
            if (options.truncate_ids)
            {
                std::ranges::copy(stream_view | view::take_until_or_throw(is_blank)
                                         | view::char_to<value_type_t<id_type>>,
                             std::back_inserter(id));
                detail::consume(stream_view | view::take_until_or_throw(is_tab));
            }
            else
            {
                std::ranges::copy(stream_view | view::take_until_or_throw(is_tab)
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
            std::ranges::next(begin(stream_view));
            detail::consume(stream_view | view::take_until_or_throw(is_tab));
        }
        std::ranges::next(begin(stream_view));

        auto seq_view{stream_view | view::take_until_or_throw(is_tab)}; // until next tab

         // Sequence
         if (is_char<'*'>(*begin(stream_view)))
             throw parse_error{"The Sequence field may not be empty for sequence files."};

         if constexpr (!detail::decays_to_ignore_v<seq_type>)
         {
             auto constexpr is_legal_alph = is_in_alphabet<seq_legal_alph_type>;
             std::ranges::copy(seq_view | std::view::transform([is_legal_alph] (char const c) // enforce legal alphabet
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

              sequence_size_after = std::ranges::size(sequence);
          }
          else
          {
              for (auto it = begin(seq_view); it != end(seq_view); ++it)
                  ++sequence_size_after;
          }
          std::ranges::next(begin(stream_view)); // skip tab

          /* Qualities */
          if (is_char<'*'>(*begin(stream_view)))
          {
              std::ranges::next(begin(stream_view)); // skip *
          }
          else
          {
             auto qual{stream_view | view::take_until_or_throw(is_tab || is_cntrl)};
             auto qual_view{qual | view::take_exactly_or_throw(sequence_size_after - sequence_size_before)};

             if constexpr (seq_qual_combined)
             {
                  // seq_qual field implies that they are the same variable
                  assert(std::addressof(sequence) == std::addressof(qualities));
                  std::ranges::copy(qual_view | view::char_to<typename value_type_t<qual_type>::quality_alphabet_type>,
                               std::ranges::begin(qualities) + sequence_size_before);
             }
             else if constexpr (!detail::decays_to_ignore_v<qual_type>)
             {
                  std::ranges::copy(qual_view | view::char_to<value_type_t<qual_type>>,
                               std::back_inserter(qualities));
             }
             else
             {
                  detail::consume(qual_view);
             }

             if (!(is_tab(*begin(qual)) | is_cntrl(*begin(qual))))
                throw unexpected_end_of_input{"Quality length surpasses sequence length."};
          }

          // consume the remaining characters (optional tags)
          detail::consume(stream_view | view::take_until_or_throw(is_cntrl));
          std::ranges::next(begin(stream_view)); //consume newline

          // make sure "buffer at end" implies "stream at end"
          if ((std::istreambuf_iterator<char>{stream} == std::istreambuf_iterator<char>{}) &&
              (!stream.eof()))
          {
              stream.get(); // triggers error in stream and sets eof
          }
    }

    //!\copydoc SequenceFileOutputFormat::write
    template <typename stream_type,     // constraints checked by file
              typename seq_type,        // other constraints checked inside function
              typename id_type,
              typename qual_type>
    void write(stream_type                        & stream,
               sequence_file_output_options const & options,
               seq_type                           && sequence,
               id_type                            && id,
               qual_type                          && qualities)
    {
        if constexpr (!(detail::decays_to_ignore_v<seq_type>))
        {
            static_assert(std::ranges::ForwardRange<seq_type> && Alphabet<value_type_t<seq_type>>,
                          "The sequence must model std::ranges::ForwardRange and its value type must model "
                          "seqan3::Alphabet.");
        }
        std::ranges::ostreambuf_iterator stream_it{stream};
        // ID
        if constexpr (detail::decays_to_ignore_v<id_type>)
            stream_it = '*';
        else if (std::ranges::empty(id)) //[[unlikely]]
            stream_it = '*';
        else
            std::ranges::copy(id, stream_it);

        stream << "\t0\t*\t0\t0\t*\t*\t0\t0\t";

        // Sequence
        if constexpr (detail::decays_to_ignore_v<seq_type>)
            stream_it = '*';
        else if (std::ranges::empty(sequence)) //[[unlikely]]
            stream_it = '*';
        else
            std::ranges::copy(sequence | view::to_char, stream_it);
        stream_it = '\t';

        // Quality line
        if constexpr (detail::decays_to_ignore_v<qual_type>)
            stream_it = '*';
        else if (std::ranges::empty(qualities))
            stream_it = '*';
        else
            std::ranges::copy(qualities | view::to_char, stream_it);

        detail::write_eol(stream_it, options.add_carriage_return);
    }
};

} // namespace seqan3
