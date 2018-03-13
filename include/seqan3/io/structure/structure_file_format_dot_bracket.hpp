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
 * \brief Provides the seqan3::structure_file_format_dot_bracket class.
 * \author Jörg Winkler <j.winkler AT fu-berlin.de>
 */

#pragma once

#include <charconv>
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
#include <range/v3/view/transform.hpp>
#include <range/v3/view/take_while.hpp>

#include <seqan3/alphabet/nucleotide/rna5.hpp>
#include <seqan3/core/metafunction/range.hpp>
#include <seqan3/io/detail/ignore_output_iterator.hpp>
#include <seqan3/io/detail/output_iterator_conversion_adaptor.hpp>
#include <seqan3/io/stream/parse_condition.hpp>
#include <seqan3/io/structure/structure_file_in_format_concept.hpp>
#include <seqan3/io/structure/structure_file_in.hpp>
#include <seqan3/range/view/char_to.hpp>
#include <seqan3/range/view/to_char.hpp>
#include <seqan3/std/concept/range.hpp>
#include <seqan3/std/view/subrange.hpp>
#include <seqan3/std/view/transform.hpp>

namespace seqan3
{

class structure_file_format_dot_bracket
{
public:
    /*!\name Constructors, destructor and assignment
     * \brief Rule of five explicitly defaulted.
     * \{
     */
    structure_file_format_dot_bracket() = default;
    structure_file_format_dot_bracket(structure_file_format_dot_bracket const &) = delete;
    structure_file_format_dot_bracket & operator=(structure_file_format_dot_bracket const &) = delete;
    structure_file_format_dot_bracket(structure_file_format_dot_bracket &&) = default;
    structure_file_format_dot_bracket & operator=(structure_file_format_dot_bracket &&) = default;
    //!\}

    //!\brief The valid file extensions for this format; note that you can modify this value.
    static inline std::vector<std::string> file_extensions
    {
        { "dbn" }
    };

    //!\copydoc sequence_file_in_format_concept::read
    template <typename stream_type,     // constraints checked by file
              typename seq_legal_alph_type,
              typename structure_legal_alph_type,
              typename seq_type,        // other constraints checked inside function
              typename id_type,
              typename bpp_type,
              typename structure_type,
              typename structured_seq_type,
              typename energy_type,
              typename react_type,
              typename comment_type,
              typename offset_type>
    void read(stream_type & stream,
              structure_file_in_options<seq_legal_alph_type, structure_legal_alph_type> const & options,
              seq_type & seq,
              id_type & id,
              bpp_type & bpp,
              structure_type & structure,
              structured_seq_type & structured_seq,
              energy_type & energy,
              react_type & react,
              react_type & react_error,
              comment_type & comment,
              offset_type & offset)
    {
        static_assert(detail::decays_to_ignore_v<seq_type> || detail::decays_to_ignore_v<structured_seq_type>,
                      "Either the sequence field, or the structured sequence field need to be set to std::ignore.");
        static_assert(detail::decays_to_ignore_v<structure_type> || detail::decays_to_ignore_v<structured_seq_type>,
                      "Either the structure field, or the structured sequence field need to be set to std::ignore.");

        auto stream_view = ranges::iterator_range{std::istreambuf_iterator<char>{stream},
                                                  std::istreambuf_iterator<char>{}};
        auto const is_eol = is_char<'\n'>{} || is_char<'\r'>{};

        // READ ID
        auto const is_id = is_char<'>'>{} || is_char<';'>{};
        if (!is_id(*ranges::begin(stream_view)))
            throw parse_error{std::string{"Expected to be on beginning of ID, but "} + is_id.msg.string() +
                              " evaluated to false on " + detail::make_printable(*ranges::begin(stream_view))};
        if (options.truncate_ids)
        {
            ranges::copy(stream_view | ranges::view::drop_while(is_id || is_blank)        // skip leading >
                                     | ranges::view::take_while(!(is_cntrl || is_blank)), // read ID until delimiter…
                         detail::make_conversion_output_iterator(id));                    // … ^A is old delimiter
            // consume rest of line
            ranges::copy(stream_view | ranges::view::take_while(!is_eol),
                         detail::make_conversion_output_iterator(std::ignore));
        }
        else
        {
            ranges::copy(stream_view | ranges::view::drop_while(is_id || is_blank)        // skip leading >
                                     | ranges::view::take_while(!is_eol),                 // read line
                         detail::make_conversion_output_iterator(id));
        }
        // skip newline
        auto it = ranges::begin(stream_view);
        if (*it == '\r')
            ++it;
        if (*it == '\n')
            ++it;

        // READ SEQUENCE
        if constexpr (!detail::decays_to_ignore_v<seq_type>)
        {
            is_in_alphabet<seq_legal_alph_type> const is_legal_alph;
            ranges::copy(stream_view | ranges::view::take_while(!is_eol)               // until end of line
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
                                     | view::char_to<value_type_t<seq_type>>,     // convert to actual target alphabet
                         detail::make_conversion_output_iterator(seq));
//            std::cout << "read sequence: " << (sequence | view::to_char ) << std::endl;
        }
        else
        {
            auto seq_view = stream_view | ranges::view::take_while(!is_eol);
            for (auto iter = ranges::begin(seq_view); iter != ranges::end(seq_view); ++iter) {}
//            std::cout << "skip sequence" << std::endl;
        }
        // skip newline
        it = ranges::begin(stream_view);
        if (*it == '\r')
            ++it;
        if (*it == '\n')
            ++it;

        // READ STRUCTURE
        if constexpr (!detail::decays_to_ignore_v<structure_type>)
        {
            is_in_alphabet<structure_legal_alph_type> const is_legal_structure;
            ranges::copy(stream_view | ranges::view::take_while(!is_space) // until whitespace
                                     | ranges::view::transform([is_legal_structure](char const c)
                                       {
                                           if (!is_legal_structure(c))
                                           {
                                               throw parse_error{
                                                   std::string{"Encountered an unexpected letter: "} +
                                                   is_legal_structure.msg.string() +
                                                   " evaluated to false on " + detail::make_printable(c)};
                                           }
                                           return c;
                                       })                                        // enforce legal alphabet
                                     | view::char_to<value_type_t<structure_type>>,    // convert to actual target alphabet
                         detail::make_conversion_output_iterator(structure));
        }
        else
        {
            auto seq_view = stream_view | ranges::view::take_while(!is_space);     // until whitespace
            for (auto iter = ranges::begin(seq_view); iter != ranges::end(seq_view); ++iter) {}
        }

        // READ ENERGY
        if constexpr (!detail::decays_to_ignore_v<energy_type>)
        {
            std::string e_str = stream_view | ranges::view::take_while(!is_eol)
                                            | ranges::view::remove_if(is_space || is_char<'('>{} || is_char<')'>{});

            // std::from_chars_result res = std::from_chars(e_str.data(), e_str.data() + e_str.size(), energy);
            // does not work (why?)

            if (!e_str.empty())
            {
                size_t num_processed = 0ul;
                energy = std::stod(e_str, &num_processed);
                if (num_processed != e_str.size())
                {
                    throw parse_error{std::string{"Failed to parse energy value '"} + e_str + "'."};
                }
            }
        }
        else
        {
            auto seq_view = stream_view | ranges::view::take_while(!is_eol);     // until newline
            for (auto iter = ranges::begin(seq_view); iter != ranges::end(seq_view); ++iter) {}
        }
        // skip newline
        it = ranges::begin(stream_view);
        if (*it == '\r')
            ++it;
        if (*it == '\n')
            ++it;

        // make sure "buffer at end" implies "stream at end"
        if ((std::istreambuf_iterator<char>{stream} == std::istreambuf_iterator<char>{}) &&
            (!stream.eof()))
        {
            stream.get(); // triggers error in stream and sets eof
        }
    }

//    template<typename stream_type,     // constraints checked by file
//             typename options_type,    // given by file
//             typename seq_type,        // other constraints checked inside function
//             typename id_type>
//    void write([[maybe_unused]] stream_type &stream,
//               [[maybe_unused]] options_type const &options,
//               [[maybe_unused]] seq_type &&seq,
//               [[maybe_unused]] id_type &&id)
//    {
//        //TODO
//    }
};

/** implementations **/


} // namespace seqan3

