// -----------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2020, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2020, Knut Reinert & MPI für molekulare Genetik
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
#include <seqan3/range/views/repeat_n.hpp>
#include <seqan3/range/views/to_char.hpp>
#include <seqan3/range/views/take_line.hpp>
#include <seqan3/range/views/take_until.hpp>
#include <seqan3/std/algorithm>
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
                std::ranges::copy(idbuffer | views::char_to<std::ranges::range_value_t<id_type>>, std::ranges::back_inserter(id));
                do
                {
                    std::ranges::copy(stream_view | views::take_until_or_throw(is_char<'S'>)
                                                  | views::char_to<std::ranges::range_value_t<id_type>>,
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
                                                  | views::char_to<std::ranges::range_value_t<id_type>>,
                                 std::ranges::back_inserter(id));
                }
                else
                {
                    std::ranges::copy(stream_view | views::take_until_or_throw(is_char<';'>)
                                                  | views::char_to<std::ranges::range_value_t<id_type>>,
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
                                  | views::char_to<std::ranges::range_value_t<seq_type>>,         // convert to actual target alphabet
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
        seqan3::detail::fast_ostreambuf_iterator stream_it{*stream.rdbuf()};

        [[maybe_unused]] size_t sequence_size = 0;

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
                stream_it.write_range(id);
            }
            else
            {
                stream_it.write_range(std::string_view{"ID "});
                stream_it.write_range(id);
                stream_it.write_range(std::string_view{"; "});
                stream_it.write_number(sequence_size);
                stream_it.write_range(std::string_view{" BP.\n"});
            }
        }

        // Sequence
        if constexpr (detail::decays_to_ignore_v<seq_type>) // sequence
        {
            throw std::logic_error{"The SEQ field may not be set to ignore when writing embl files."};
        }
        else
        {
            if (std::ranges::empty(sequence)) //[[unlikely]]
                throw std::runtime_error{"The SEQ field may not be empty when writing embl files."};

            // write beginning of sequence record
            stream_it.write_range(std::string_view{"SQ Sequence "});
            stream_it.write_number(sequence_size);
            stream_it.write_range(std::string_view{" BP;\n"});

            // write sequence in chunks of 60 bp's with a space after 10 bp's
            auto char_sequence = sequence | views::to_char;
            auto it = std::ranges::begin(char_sequence);
            size_t written_chars{0};
            uint8_t chunk_size{10u};

            while (it != std::ranges::end(char_sequence))
            {
                auto current_end = it;
                size_t steps = std::ranges::advance(current_end, chunk_size, std::ranges::end(char_sequence));

                using subrange_t = std::ranges::subrange<decltype(it), decltype(it), std::ranges::subrange_kind::sized>;
                it = stream_it.write_range(subrange_t{it, current_end, chunk_size - steps});
                stream_it = ' ';
                written_chars += chunk_size;

                if (written_chars % 60 == 0)
                {
                    stream_it.write_number(written_chars);
                    stream_it.write_end_of_line(options.add_carriage_return);
                }
            }

            // fill last line
            auto characters_in_last_line = sequence_size % 60;
            auto number_of_padding_needed = 65 - characters_in_last_line - characters_in_last_line / chunk_size;
            stream_it.write_range(views::repeat_n(' ', number_of_padding_needed));
            stream_it.write_number(sequence_size);
            stream_it.write_end_of_line(options.add_carriage_return);

            // write end-of-record-symbol
            stream_it.write_range(std::string_view{"//"});
            stream_it.write_end_of_line(options.add_carriage_return);
        }
    }
};

} // namespace seqan
