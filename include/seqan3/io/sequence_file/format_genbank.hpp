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
 * \brief Provides the seqan3::sequence_file_format_genbank class.
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
/*!\brief       The GenBank format.
 * \implements  sequence_file_input_format_concept
 * \implements  sequence_file_output_format_concept
 * \ingroup     sequence
 *
 * \details
 *
 * ### Introduction
 *
 * genbank is the format used in GenBank sequence database. See
 * (http://quma.cdb.riken.jp/help/gbHelp.html) for a an in-depth description of the format.
 *
 * ### Fields
 *
 * The genbank format provides the fields seqan3::field::SEQ and seqan3::field::ID. Both fields are required when
 * writing.
 *
 * ### Implementation notes
 *
 * There is no truncate_ids option while reading because the GenBank format has no (FASTA-like) identifier. Instead,
 * there is the option t"complete_header" to indicate whether the whole header is to be read (complete_header=true) or
 * only the "LOCUS" information should be stored.
 *
 * Passed qualities to either the read or write function are ignored.
 *
 */
class sequence_file_format_genbank
{
public:
    /*!\name Constructors, destructor and assignment
     * \{
     */
     //!\brief Default constructor is explicitly deleted, you need to provide a stream or file name.
    sequence_file_format_genbank() = default;
    //!\brief Copy construction is explicitly deleted, because you can't have multiple access to the genbank file.
    sequence_file_format_genbank(sequence_file_format_genbank const &) = delete;
    //!\brief Copy assignment is explicitly deleted, because you can't have multiple access to the genbank file.
    sequence_file_format_genbank & operator=(sequence_file_format_genbank const &) = delete;
    //!\brief Move construction is defaulted.
    sequence_file_format_genbank(sequence_file_format_genbank &&) = default;
    //!\brief Move assignment is defaulted.
    sequence_file_format_genbank & operator=(sequence_file_format_genbank &&) = default;
    //!\}

    //!\brief The valid file extensions for this format; note that you can modify this value.
    static inline std::vector<std::string> file_extensions
    {
        { "genbank" },
        { "gb" },
        { "gbk" },
    };

    //!\copydoc sequence_file_input_format_concept::read
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
        auto stream_view = view::subrange<decltype(std::istreambuf_iterator<char>{stream}),
                                          decltype(std::istreambuf_iterator<char>{})>
                           {std::istreambuf_iterator<char>{stream},
                            std::istreambuf_iterator<char>{}};
        auto stream_it = ranges::begin(stream_view);

        std::string identifier;
        ranges::copy(stream_view | view::take_until_or_throw(is_cntrl || is_blank),
                     std::back_inserter(identifier));

        if (identifier != "LOCUS")
            throw parse_error{"An entry has to start with the code word LOCUS."};

        //ID
        if constexpr (!detail::decays_to_ignore_v<id_type>)
        {
            if (options.complete_header)
            {
                while (identifier != "ORIGIN")
                {
                        ranges::copy(identifier  | view::char_to<value_type_t<id_type>>,
                                     std::back_inserter(id));
                        ranges::copy(stream_view | view::take_until_or_throw(is_char<'\n'>)
                                                 | view::char_to<value_type_t<id_type>>,
                                     std::back_inserter(id));
                        ranges::copy(std::string_view{"\n"} | view::char_to<value_type_t<id_type>>,
                                     std::back_inserter(id));

                    ++stream_it;
                    identifier.clear();
                    ranges::copy(stream_view | view::take_until_or_throw(is_cntrl || is_blank),
                                 std::back_inserter(identifier));
                }
            }
            else
            {
                detail::consume(stream_view | ranges::view::take_while(is_blank));
                // read id
                ranges::copy(stream_view | view::take_until_or_throw(is_cntrl || is_blank)
                                         | view::char_to<value_type_t<id_type>>,
                             std::back_inserter(id));
                detail::consume(stream_view | view::take_line_or_throw);
                identifier.clear();
                // Update identifier to jump to sequence (see below)
                ranges::copy(stream_view | view::take_until_or_throw(is_cntrl || is_blank),
                             std::back_inserter(identifier));
            }
        }

        // Jump to sequence
        while (identifier != "ORIGIN") //True, if not complete header is read into id
        {
            detail::consume(stream_view | view::take_until(is_cntrl));
            ++stream_it;
            identifier.clear();
            ranges::copy(stream_view | view::take_until_or_throw(is_cntrl || is_blank),
                         std::back_inserter(identifier));
        }

        // Sequence
        detail::consume(stream_view | view::take_line_or_throw);
        auto constexpr is_end = is_char<'/'> ;
        if constexpr (!detail::decays_to_ignore_v<seq_type>)
        {
            auto seq_view = stream_view | ranges::view::remove_if(is_space || is_digit)
            // ignore whitespace and numbers
                                        | view::take_until_or_throw(is_end);   // until //

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
                                  | view::char_to<value_type_t<seq_type>>,         // convert to actual target alphabet
                         std::back_inserter(sequence));
        }
        else
        {
            detail::consume(stream_view | view::take_until(is_end));
        }
        //Jump over //
        detail::consume(stream_view | view::take_until(is_cntrl));
        ++stream_it;


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
               qual_type                      && SEQAN3_DOXYGEN_ONLY(qualities))
    {

        ranges::ostreambuf_iterator stream_it{stream};
        size_t sequence_size = 0;
        if constexpr (!detail::decays_to_ignore_v<seq_type>)
            sequence_size = ranges::size(sequence);

        // ID
        if constexpr (detail::decays_to_ignore_v<id_type>)
        {
            throw std::logic_error{"The ID field may not be set to ignore when writing genbank files."};
        }
        else
        {
            if (ranges::empty(id)) //[[unlikely]]
                throw std::runtime_error{"The ID field may not be empty when writing genbank files."};

            if (options.complete_header)
            {
                ranges::copy(id, stream_it);
            }
            else
            {
                ranges::copy(std::string_view{"LOCUS       "}, stream_it);
                ranges::copy(id, stream_it);
                ranges::copy(std::string_view{"                 "}, stream_it);
                ranges::copy(std::to_string(sequence_size), stream_it);
                ranges::copy(std::string_view{" bp\n"}, stream_it);
            }
        }

        // Sequence
        if constexpr (detail::decays_to_ignore_v<seq_type>) // sequence
        {
            throw std::logic_error{"The SEQ and SEQ_QUAL fields may not both be set to ignore when writing genbank files."};
        }
        else
        {
            if (ranges::empty(sequence)) //[[unlikely]]
                throw std::runtime_error{"The SEQ field may not be empty when writing genbank files."};

            ranges::copy(std::string_view{"ORIGIN"}, stream_it);
            stream_it = '\n';
            auto seq = sequence | ranges::view::chunk(60);
            unsigned int i = 0;
            size_t bp = 1;
            while (bp < sequence_size)
            {
                for(size_t j = std::to_string(bp).length(); j < 9; j++)
                    stream_it = ' ';
                ranges::copy(std::to_string(bp), stream_it);
                stream_it = ' ';
                ranges::copy(seq[i] | view::to_char
                                    | ranges::view::chunk(10)
                                    | ranges::view::join(' '), stream_it);
                bp += 60;
                ++i;
                stream_it = '\n';
            }
            ranges::copy(std::string_view{"//"}, stream_it);
            stream_it = '\n';
        }
    }

};

} // namespace seqan3
