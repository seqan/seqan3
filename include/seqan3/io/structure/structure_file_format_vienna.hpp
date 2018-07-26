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
 * \brief Provides the seqan3::structure_file_format_vienna class.
 * \author Jörg Winkler <j.winkler AT fu-berlin.de>
 */

#pragma once

#include <iterator>
#include <stack>
#include <string>
#include <string_view>
#include <type_traits>
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
#include <seqan3/io/detail/misc.hpp>
#include <seqan3/io/stream/parse_condition.hpp>
#include <seqan3/io/structure/detail.hpp>
#include <seqan3/io/structure/structure_file_in_options.hpp>
#include <seqan3/io/structure/structure_file_out_options.hpp>
#include <seqan3/range/detail/misc.hpp>
#include <seqan3/range/view/char_to.hpp>
#include <seqan3/range/view/to_char.hpp>
#include <seqan3/range/view/take_line.hpp>
#include <seqan3/std/concept/range.hpp>
#include <seqan3/std/view/subrange.hpp>
#include <seqan3/std/view/transform.hpp>

namespace seqan3
{
/*!\brief       The Vienna format (dot bracket notation) for RNA sequences with secondary structure.
 * \implements  seqan3::structure_file_in_format_concept
 * \implements  seqan3::structure_file_out_format_concept
 * \ingroup     structure
 *
 * \details
 *
 * ### Introduction
 *
 * Dot Bracket or Vienna Notation is widely used for secondary structure annotation. Is is a very simple format,
 * containing one or more sequences. Each sequence must appear as a single line in the file.
 * A sequence may be preceded by a special line starting with the '>' character followed by a sequence name
 * (like FastA). After each sequence line there is usually a line containing
 * secondary structure, using brackets to denote interacting nucleotides or amino acids, and dots for unpaired sites.
 * The length of the struture must equal the length of the sequence.
 * Optionally, the structure may be followed by a space character and the minimum free energy value enclosed
 * in parentheses (). Note that there cannot be energy without structure.
 *
 * The Vienna format is the output format of _RNAfold_. Furthermore, it is designed to be compatible with the
 * input format of the ViennaRNA package (if structure and energy are omitted).
 * See https://www.tbi.univie.ac.at/RNA/tutorial/#sec2_7 for details.
 *
 * ### Fields
 *
 * The Vienna format provides the fields seqan3::field::SEQ, seqan3::field::ID, seqan3::field::STRUCTURE,
 * seqan3::field::STRUCTURED_SEQ and seqan3::field::ENERGY.\n
 * If you select seqan3::field::STRUCTURED_SEQ you must not select seqan3::field::SEQ or seqan3::field::STRUCTURE.
 * Either the field seqan3::field::SEQ or the field seqan3::field::STRUCTURED_SEQ is required when writing.
 *
 * ### Implementation notes
 *
 * In the unlikely case that you want to read a Vienna file, which does not contain structure lines, you have to
 * assign false to `seqan3::structure_file_in_options::file_has_structure`. \n\n
 * When reading the ID-line the identifier (`>`) and any blank characters before the actual ID are
 * stripped. Each field is read/written as a single line (except ENERGY, which goes right after the structure).
 * Numbers and spaces within the sequence are simply ignored, but not within the structure.
 */
class structure_file_format_vienna
{
public:
    /*!\name Constructors, destructor and assignment
     * \brief Rule of five explicitly defaulted.
     * \{
     */
    structure_file_format_vienna() = default;
    structure_file_format_vienna(structure_file_format_vienna const &) = delete;
    structure_file_format_vienna & operator=(structure_file_format_vienna const &) = delete;
    structure_file_format_vienna(structure_file_format_vienna &&) = default;
    structure_file_format_vienna & operator=(structure_file_format_vienna &&) = default;
    //!\}

    //!\brief The valid file extensions for this format; note that you can modify this value.
    static inline std::vector<std::string> file_extensions
    {
        { "dbn" }
    };

    //!\copydoc seqan3::structure_file_in_format_concept::read
    template <typename stream_type,     // constraints checked by file
              typename seq_legal_alph_type,
              bool     structured_seq_combined,
              typename seq_type,        // other constraints checked inside function
              typename id_type,
              typename bpp_type,
              typename structure_type,
              typename energy_type,
              typename react_type,
              typename comment_type,
              typename offset_type>
    void read(stream_type & stream,
              structure_file_in_options<seq_legal_alph_type, structured_seq_combined> const & options,
              seq_type & seq,
              id_type & id,
              bpp_type & bpp,
              structure_type & structure,
              energy_type & energy,
              react_type & SEQAN3_DOXYGEN_ONLY(react),
              react_type & SEQAN3_DOXYGEN_ONLY(react_err),
              comment_type & SEQAN3_DOXYGEN_ONLY(comment),
              offset_type & SEQAN3_DOXYGEN_ONLY(offset))
    {
        auto stream_view = view::subrange<decltype(std::istreambuf_iterator<char>{stream}),
                                          decltype(std::istreambuf_iterator<char>{})>
                           { std::istreambuf_iterator<char>{stream}, std::istreambuf_iterator<char>{} };
        // READ ID (if present)
        auto const is_id = is_char<'>'>{};
        if (is_id(*ranges::begin(stream_view)))
        {
            if constexpr (!detail::decays_to_ignore_v<id_type>)
            {
                if (options.truncate_ids)
                {
                    ranges::copy(stream_view | ranges::view::drop_while(is_id || is_blank)        // skip leading >
                                 | ranges::view::take_while(!(is_cntrl || is_blank)), // read ID until delimiter
                                 detail::make_conversion_output_iterator(id));                    // … ^A is old delimiter
                    detail::consume(stream_view | view::take_line_or_throw);
                }
                else
                {
                    ranges::copy(stream_view | ranges::view::drop_while(is_id || is_blank)        // skip leading >
                                 | view::take_line_or_throw,                          // read line
                                 detail::make_conversion_output_iterator(id));
                }
            }
            else
            {
                detail::consume(stream_view | view::take_line_or_throw);
            }
        }

        // READ SEQUENCE
        if constexpr (!detail::decays_to_ignore_v<seq_type>)
        {
            is_in_alphabet<seq_legal_alph_type> const is_legal_seq;
            ranges::copy(stream_view | view::take_line_or_throw                      // until end of line
                                     | ranges::view::remove_if(is_space || is_digit) // ignore whitespace and numbers
                                     | ranges::view::transform([is_legal_seq](char const c)
                                       {
                                           if (!is_legal_seq(c))                    // enforce legal alphabet
                                           {
                                               throw parse_error{std::string{"Encountered an unexpected letter: "} +
                                                                 is_legal_seq.msg.string() +
                                                                 " evaluated to false on " +
                                                                 detail::make_printable(c)};
                                           }
                                         return c;
                                       })
                                     | view::char_to<value_type_t<seq_type>>, // convert to actual target alphabet
                         detail::make_conversion_output_iterator(seq));
        }
        else
        {
            detail::consume(stream_view | view::take_line);
            detail::consume(stream_view | ranges::view::take_while(is_space));
        }

        if (!options.file_has_structure)
        {
            if ((std::istreambuf_iterator<char>{stream} == std::istreambuf_iterator<char>{}) && (!stream.eof()))
                stream.get();

            return;
        }

        // READ STRUCTURE
        if constexpr (!detail::decays_to_ignore_v<structure_type>)
        {
            if constexpr (structured_seq_combined)
            {
                assert(std::addressof(seq) == std::addressof(structure));
                using alph_type = typename value_type_t<structure_type>::structure_alphabet_type;
                ranges::copy(read_structure<alph_type>(stream_view), ranges::begin(structure));

                if constexpr (!detail::decays_to_ignore_v<bpp_type>)
                    detail::bpp_from_rna_structure<alph_type>(bpp, structure);
            }
            else
            {
                using alph_type = value_type_t<structure_type>;
                ranges::copy(read_structure<alph_type>(stream_view),
                             detail::make_conversion_output_iterator(structure));

                if constexpr (!detail::decays_to_ignore_v<bpp_type>)
                    detail::bpp_from_rna_structure<alph_type>(bpp, structure);
            }
        }
        else if constexpr (!detail::decays_to_ignore_v<bpp_type>)
        {
            detail::bpp_from_rna_structure<wuss51>(bpp, read_structure<wuss51>(stream_view));
        }
        else
        {
            detail::consume(stream_view | ranges::view::take_while(!is_space)); // until whitespace
        }

        // READ ENERGY (if present)
        if constexpr (!detail::decays_to_ignore_v<energy_type>)
        {
            std::string e_str = stream_view | view::take_line
                                            | ranges::view::remove_if(is_space || is_char<'('>{} || is_char<')'>{});
            if (!e_str.empty())
            {
                size_t num_processed;
                energy = std::stod(e_str, &num_processed);
                if (num_processed != e_str.size()) // [[unlikely]]
                {
                    throw parse_error{std::string{"Failed to parse energy value '"} + e_str + "'."};
                }
            }
        }
        else
        {
            detail::consume(stream_view | view::take_line);
            detail::consume(stream_view | ranges::view::take_while(is_space));
        }

        // make sure "buffer at end" implies "stream at end"
        if ((std::istreambuf_iterator<char>{stream} == std::istreambuf_iterator<char>{}) && (!stream.eof()))
        {
            stream.get(); // triggers error in stream and sets eof
        }
    }

    //!\copydoc seqan3::structure_file_out_format_concept::write
    template <typename stream_type,     // constraints checked by file
              typename seq_type,        // other constraints checked inside function
              typename id_type,
              typename bpp_type,
              typename structure_type,
              typename energy_type,
              typename react_type,
              typename comment_type,
              typename offset_type>
    void write(stream_type & stream,
               structure_file_out_options const & options,
               seq_type && seq,
               id_type && id,
               bpp_type && bpp,
               structure_type && structure,
               energy_type && energy,
               react_type && SEQAN3_DOXYGEN_ONLY(react),
               react_type && SEQAN3_DOXYGEN_ONLY(react_err),
               comment_type && SEQAN3_DOXYGEN_ONLY(comment),
               offset_type && SEQAN3_DOXYGEN_ONLY(offset))
    {
        ranges::ostreambuf_iterator stream_it{stream};

        // WRITE ID (optional)
        if constexpr (!detail::decays_to_ignore_v<id_type>)
        {
            if (!ranges::empty(id))
            {
                stream_it = '>';
                stream_it = ' ';
                ranges::copy(id, stream_it);
                detail::write_eol(stream_it, options.add_carriage_return);
            }
        }

        // WRITE SEQUENCE
        if constexpr (!detail::decays_to_ignore_v<seq_type>)
        {
            if (ranges::empty(seq)) //[[unlikely]]
                throw std::runtime_error{"The SEQ field may not be empty when writing Vienna files."};

            ranges::copy(seq | view::to_char, stream_it);
            detail::write_eol(stream_it, options.add_carriage_return);
        }
        else
        {
            throw std::logic_error{"The SEQ and STRUCTURED_SEQ fields may not both be set to ignore "
                                   "when writing Vienna files."};
        }

        // WRITE STRUCTURE (optional)
        if constexpr (!detail::decays_to_ignore_v<structure_type>)
        {
            if (!ranges::empty(structure))
                ranges::copy(structure | view::to_char, stream_it);

            // WRITE ENERGY (optional)
            if constexpr (!detail::decays_to_ignore_v<energy_type>)
            {
                if (energy)
                {
                    stream_it = ' ';
                    stream_it = '(';
                    ranges::copy(std::to_string(energy), stream_it);
                    stream_it = ')';
                }
            }
            detail::write_eol(stream_it, options.add_carriage_return);
        }
        else if constexpr (!detail::decays_to_ignore_v<energy_type>)
        {
            throw std::logic_error{"The ENERGY field cannot be written to a Vienna file without providing STRUCTURE."};
        }
    }

private:
    /*!
     * \brief Extract the structure string from the given stream.
     * \tparam alph_type        The alphabet type the structure is converted to.
     * \tparam stream_view_type The type of the input stream.
     * \param stream_view       The input stream to be read.
     * \return                  A ranges::view containing the structure annotation string.
     */
    template <typename alph_type, typename stream_view_type>
    auto read_structure(stream_view_type & stream_view)
    {
        is_in_alphabet<alph_type> const is_legal_structure;
        return stream_view | ranges::view::take_while(!is_space) // until whitespace
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
                             })                                  // enforce legal alphabet
                           | view::char_to<alph_type>;           // convert to actual target alphabet
    }
};

} // namespace seqan3
