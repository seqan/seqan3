// -----------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2019, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2019, Knut Reinert & MPI für molekulare Genetik
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at: https://github.com/seqan/seqan3/blob/master/LICENSE.md
// -----------------------------------------------------------------------------------------------------

/*!\file
 * brief Provides the seqan3::format_embl.
 * \author Mitra Darvish <mitra.darvish AT fu-berlin.de>
 */

#pragma once

#include <iterator>
#include <string>
#include <string_view>
#include <vector>

#include <range/v3/view/chunk.hpp>

#include <seqan3/alphabet/nucleotide/dna5.hpp>
#include <seqan3/core/char_operations/predicate.hpp>
#include <seqan3/core/type_traits/range.hpp>
#include <seqan3/io/detail/misc.hpp>
#include <seqan3/io/sequence_file/input_format_concept.hpp>
#include <seqan3/io/sequence_file/input_options.hpp>
#include <seqan3/io/sequence_file/output_format_concept.hpp>
#include <seqan3/io/sequence_file/output_options.hpp>
#include <seqan3/io/stream/iterator.hpp>
#include <seqan3/range/detail/misc.hpp>
#include <seqan3/range/views/char_to.hpp>
#include <seqan3/range/views/istreambuf.hpp>
#include <seqan3/range/views/join.hpp>
#include <seqan3/range/views/repeat_n.hpp>
#include <seqan3/range/views/to_char.hpp>
#include <seqan3/range/views/take_exactly.hpp>
#include <seqan3/range/views/take_line.hpp>
#include <seqan3/range/views/take_until.hpp>
#include <seqan3/range/views/to.hpp>
#include <seqan3/std/algorithm>
#include <seqan3/std/charconv>
#include <seqan3/std/ranges>

namespace seqan3
{
/*!\brief       The EMBL format.
 * \implements  seqan3::sequence_file_input_format
 * \implements  seqan3::sequence_file_output_format
 * \ingroup     sequence
 *
 * \details
 *
 * ### Introduction
 *
 * embl is the format used in ENA sequence records. See the
 * (ftp://ftp.ebi.ac.uk/pub/databases/embl/doc/usrman.txt) for a an in-depth description of the format.
 *
 * ### fields_specialisation
 *
 * The embl format provides the fields seqan3::field::seq and seqan3::field::id. Both fields are required when writing.
 *
 * ### Implementation notes
 *
 * When reading the ID-line, the ID is read until the stream encounters a ';'. Unless, the option truncate_ids is set to
 * true, then the id is read until it either sees a blank, ';' or a new line. When the option
 * "embl_genbank_complete_header" is set to true (Default: false) the whole header is read into the id.
 *
 * When writing the ID-line, the sequence length is appended.
 *
 * All other identifiers apart from 'ID' and 'SQ' are currently ignored.
 *
 * Passed qualities to either the read or write function are ignored.
 *
 */
class format_embl
{
public:
    /*!\name Constructors, destructor and assignment
     * \{
     */
    format_embl() noexcept = default; //!< Defaulted.
    format_embl(format_embl const &) noexcept = default; //!< Defaulted.
    format_embl & operator=(format_embl const &) noexcept = default; //!< Defaulted.
    format_embl(format_embl &&) noexcept = default; //!< Defaulted.
    format_embl & operator=(format_embl &&) noexcept = default; //!< Defaulted.
    ~format_embl() noexcept = default; //!< Defaulted.
    //!\}

    //!\brief The valid file extensions for this format; note that you can modify this value.
    static inline std::vector<std::string> file_extensions
    {
        { "embl" },
    };

protected:
    //!\copydoc sequence_file_input_format::read_sequence_record
    template <typename stream_type,     // constraints checked by file
              typename seq_legal_alph_type, bool seq_qual_combined,
              typename seq_type,        // other constraints checked inside function
              typename id_type,
              typename qual_type>
    void read_sequence_record(stream_type & stream,
                              sequence_file_input_options<seq_legal_alph_type, seq_qual_combined> const & options,
                              seq_type    & sequence,
                              id_type     & id,
                              qual_type   & SEQAN3_DOXYGEN_ONLY(qualities))
    {
        auto stream_view = views::istreambuf(stream);
        auto stream_it = std::ranges::begin(stream_view);

        std::string idbuffer;
        std::ranges::copy(stream_view | views::take_until_or_throw(is_cntrl || is_blank),
                          std::ranges::back_inserter(idbuffer));
        if (idbuffer != "ID")
            throw parse_error{"An entry has to start with the code word ID."};

        if constexpr (!detail::decays_to_ignore_v<id_type>)
        {
            if (options.embl_genbank_complete_header)
            {
                std::ranges::copy(idbuffer | views::char_to<value_type_t<id_type>>, std::ranges::back_inserter(id));
                do
                {
                    std::ranges::copy(stream_view | views::take_until_or_throw(is_char<'S'>)
                                                  | views::char_to<value_type_t<id_type>>,
                                 std::ranges::back_inserter(id));
                    id.push_back(*stream_it);
                    ++stream_it;
                } while (*stream_it != 'Q');
                id.pop_back(); // remove 'S' from id
                idbuffer = "SQ";
            }
            else
            {
                // ID
                detail::consume(stream_view | views::take_until(!is_blank));

                // read id
                if (options.truncate_ids)
                {
                    std::ranges::copy(stream_view | views::take_until_or_throw(is_blank || is_char<';'> || is_cntrl)
                                                  | views::char_to<value_type_t<id_type>>,
                                 std::ranges::back_inserter(id));
                }
                else
                {
                    std::ranges::copy(stream_view | views::take_until_or_throw(is_char<';'>)
                                                  | views::char_to<value_type_t<id_type>>,
                                 std::ranges::back_inserter(id));
                }
            }
        }

        // Jump to sequence
        if (idbuffer !="SQ")
        {
            do
            {
                detail::consume(stream_view | views::take_until_or_throw(is_char<'S'>));
                ++stream_it;
            } while (*stream_it != 'Q');
        }
        detail::consume(stream_view | views::take_line_or_throw); //Consume line with infos to sequence

        // Sequence
        auto constexpr is_end = is_char<'/'> ;
        if constexpr (!detail::decays_to_ignore_v<seq_type>)
        {
            auto seq_view = stream_view | std::views::filter(!(is_space || is_digit)) // ignore whitespace and numbers
                                        | views::take_until_or_throw(is_end);   // until //

            auto constexpr is_legal_alph = is_in_alphabet<seq_legal_alph_type>;
            std::ranges::copy(seq_view | std::views::transform([is_legal_alph] (char const c) // enforce legal alphabet
                                    {
                                        if (!is_legal_alph(c))
                                        {
                                            throw parse_error{std::string{"Encountered an unexpected letter: "} +
                                                              is_legal_alph.msg +
                                                              " evaluated to false on " +
                                                              detail::make_printable(c)};
                                        }
                                        return c;
                                    })
                                  | views::char_to<value_type_t<seq_type>>,         // convert to actual target alphabet
                         std::ranges::back_inserter(sequence));
        }
        else
        {
            detail::consume(stream_view | views::take_until(is_end));
        }
        //Jump over // and cntrl
        ++stream_it;
        ++stream_it;
        ++stream_it;
    }

    //!\copydoc sequence_file_output_format::write_sequence_record
    template <typename stream_type,     // constraints checked by file
              typename seq_type,        // other constraints checked inside function
              typename id_type,
              typename qual_type>
    void write_sequence_record(stream_type                          & stream,
                               sequence_file_output_options const   & options,
                               seq_type                             && sequence,
                               id_type                              && id,
                               qual_type                            && SEQAN3_DOXYGEN_ONLY(qualities))
    {
        seqan3::ostreambuf_iterator stream_it{stream};
        [[maybe_unused]] size_t sequence_size = 0;
        [[maybe_unused]] char buffer[50];
        if constexpr (!detail::decays_to_ignore_v<seq_type>)
            sequence_size = std::ranges::distance(sequence);

        // ID
        if constexpr (detail::decays_to_ignore_v<id_type>)
        {
            throw std::logic_error{"The ID field may not be set to ignore when writing embl files."};
        }
        else
        {
            if (ranges::empty(id)) //[[unlikely]]
                throw std::runtime_error{"The ID field may not be empty when writing embl files."};

            if (options.embl_genbank_complete_header)
            {
                std::ranges::copy(id, stream_it);
            }
            else
            {
                std::ranges::copy(std::string_view{"ID "}, stream_it);
                std::ranges::copy(id, stream_it);
                std::ranges::copy(std::string_view{"; "}, stream_it);
                auto res = std::to_chars(&buffer[0], &buffer[0] + sizeof(buffer), sequence_size);
                std::copy(&buffer[0], res.ptr, stream_it);
                std::ranges::copy(std::string_view{" BP.\n"}, stream_it);
            }

        }

        // Sequence
        if constexpr (detail::decays_to_ignore_v<seq_type>) // sequence
        {
            throw std::logic_error{"The SEQ field may not be set to ignore when writing embl files."};
        }
        else
        {
            if (ranges::empty(sequence)) //[[unlikely]]
                throw std::runtime_error{"The SEQ field may not be empty when writing embl files."};

            std::ranges::copy(std::string_view{"SQ Sequence "}, stream_it);
            auto res = std::to_chars(&buffer[0], &buffer[0] + sizeof(buffer), sequence_size);
            std::copy(&buffer[0], res.ptr, stream_it);
            std::ranges::copy(std::string_view{" BP;\n"}, stream_it);
            auto seqChunk = sequence | ranges::view::chunk(60);
            unsigned int i = 0;
            size_t bp = 0;
            for (auto chunk : seqChunk)
            {
                std::ranges::copy(chunk | views::to_char
                                        | ranges::view::chunk(10)
                                        | views::join(' '), stream_it);
                ++i;
                stream_it = ' ';
                bp = std::min(sequence_size, bp + 60);
                uint8_t num_blanks = 60 * i - bp;  // for sequence characters
                num_blanks += num_blanks / 10;     // additional chunk separators
                std::ranges::copy(views::repeat_n(' ', num_blanks), stream_it);
                std::ranges::copy(std::to_string(bp), stream_it);
                stream_it = '\n';
            }
            std::ranges::copy(std::string_view{"//"}, stream_it);
            stream_it = '\n';
        }
    }
};

} // namespace seqan
