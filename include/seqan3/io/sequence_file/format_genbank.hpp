// SPDX-FileCopyrightText: 2006-2025 Knut Reinert & Freie Universität Berlin
// SPDX-FileCopyrightText: 2016-2025 Knut Reinert & MPI für molekulare Genetik
// SPDX-License-Identifier: BSD-3-Clause

/*!\file
 * \brief Provides the seqan3::sequence_file_format_genbank class.
 * \author Mitra Darvish <mitra.darvish AT fu-berlin.de>
 */

#pragma once

#include <algorithm>
#include <iterator>
#include <ranges>
#include <seqan3/std/charconv>
#include <string>
#include <string_view>
#include <vector>

#include <seqan3/alphabet/concept.hpp>
#include <seqan3/alphabet/nucleotide/dna5.hpp>
#include <seqan3/alphabet/views/char_to.hpp>
#include <seqan3/alphabet/views/to_char.hpp>
#include <seqan3/core/range/detail/misc.hpp>
#include <seqan3/core/range/type_traits.hpp>
#include <seqan3/io/detail/misc.hpp>
#include <seqan3/io/sequence_file/input_format_concept.hpp>
#include <seqan3/io/sequence_file/input_options.hpp>
#include <seqan3/io/sequence_file/output_format_concept.hpp>
#include <seqan3/io/sequence_file/output_options.hpp>
#include <seqan3/io/views/detail/istreambuf_view.hpp>
#include <seqan3/io/views/detail/take_line_view.hpp>
#include <seqan3/io/views/detail/take_until_view.hpp>
#include <seqan3/utility/char_operations/predicate.hpp>
#include <seqan3/utility/detail/type_name_as_string.hpp>
#include <seqan3/utility/views/chunk.hpp>
#include <seqan3/utility/views/interleave.hpp>

namespace seqan3
{

/*!\brief The GenBank format.
 * \implements seqan3::sequence_file_input_format
 * \implements seqan3::sequence_file_output_format
 * \ingroup io_sequence_file
 *
 * \details
 *
 * ### Introduction
 *
 * genbank is the format used in the GenBank sequence database. See [this example record at NCBI]
 * (https://www.ncbi.nlm.nih.gov/Sitemap/samplerecord.html) for more details about the format.
 *
 * ### fields_specialisation
 *
 * The genbank format provides the fields seqan3::field::seq and seqan3::field::id. Both fields are required when
 * writing.
 *
 * ### Implementation notes
 *
 * There is no truncate_ids option while reading because the GenBank format has no (FASTA-like) idbuffer. Instead,
 * there is the option "complete_header" to indicate whether the whole header is to be read
 * (embl_genbank_complete_header=true) or only the "LOCUS" information should be stored.
 *
 * Qualities passed to the write function are ignored.
 *
 * \remark For a complete overview, take a look at \ref io_sequence_file
 */
class format_genbank
{
public:
    /*!\name Constructors, destructor and assignment
     * \{
     */
    format_genbank() noexcept = default;                                   //!< Defaulted.
    format_genbank(format_genbank const &) noexcept = default;             //!< Defaulted.
    format_genbank & operator=(format_genbank const &) noexcept = default; //!< Defaulted.
    format_genbank(format_genbank &&) noexcept = default;                  //!< Defaulted.
    format_genbank & operator=(format_genbank &&) noexcept = default;      //!< Defaulted.
    ~format_genbank() noexcept = default;                                  //!< Defaulted.

    //!\}

    //!\brief The valid file extensions for this format; note that you can modify this value.
    static inline std::vector<std::string> file_extensions{
        {"genbank"},
        {"gb"},
        {"gbk"},
    };

protected:
    //!\copydoc sequence_file_input_format::read_sequence_record
    template <typename stream_type, // constraints checked by file
              typename seq_legal_alph_type,
              typename stream_pos_type,
              typename seq_type, // other constraints checked inside function
              typename id_type,
              typename qual_type>
    void read_sequence_record(stream_type & stream,
                              sequence_file_input_options<seq_legal_alph_type> const & options,
                              stream_pos_type & position_buffer,
                              seq_type & sequence,
                              id_type & id,
                              qual_type & SEQAN3_DOXYGEN_ONLY(qualities))
    {
        // Store current position in buffer
        // Must happen before constructing the view.
        // With libc++, tellg invalidates the I/O buffer.
        position_buffer = stream.tellg();

        auto stream_view = detail::istreambuf(stream);
        auto stream_it = std::ranges::begin(stream_view);

        if (!(std::ranges::equal(stream_view | detail::take_until_or_throw(is_cntrl || is_blank),
                                 std::string{"LOCUS"})))
            throw parse_error{"An entry has to start with the code word LOCUS."};

        //ID
        if constexpr (!detail::decays_to_ignore_v<id_type>)
        {
            if (options.embl_genbank_complete_header)
            {
                std::ranges::copy(std::string_view{"LOCUS"}, std::back_inserter(id));

                while (!is_char<'O'>(*std::ranges::begin(stream_view)))
                {
                    std::ranges::copy(stream_view | detail::take_line_or_throw
                                          | views::char_to<std::ranges::range_value_t<id_type>>,
                                      std::back_inserter(id));
                    id.push_back('\n');
                }
            }
            else
            {
                detail::consume(stream_view | detail::take_until(!is_blank));

                auto read_id_until = [&stream_view, &id](auto predicate)
                {
                    std::ranges::copy(stream_view | detail::take_until_or_throw(predicate)
                                          | views::char_to<std::ranges::range_value_t<id_type>>,
                                      std::back_inserter(id));
                };

                if (options.truncate_ids)
                    read_id_until(is_space);
                else
                    read_id_until(is_cntrl);

                detail::consume(stream_view | detail::take_line_or_throw);
            }
        }

        // Jump to sequence
        while (!(is_char<'O'>(*std::ranges::begin(stream_view)) || options.embl_genbank_complete_header))
            detail::consume(stream_view | detail::take_line_or_throw);

        // Sequence
        detail::consume(stream_view | detail::take_line_or_throw); // consume "ORIGIN"
        constexpr auto is_end = is_char<'/'>;
        if constexpr (!detail::decays_to_ignore_v<seq_type>)
        {
            constexpr auto is_legal_alph = char_is_valid_for<seq_legal_alph_type>;
            std::ranges::copy(
                stream_view | std::views::filter(!(is_space || is_digit)) // ignore whitespace and numbers
                    | detail::take_until_or_throw(is_end)                 // until //
                    | std::views::transform(
                        [is_legal_alph](char const c) // enforce legal alphabet
                        {
                            if (!is_legal_alph(c))
                            {
                                throw parse_error{std::string{"Encountered an unexpected letter: "}
                                                  + "char_is_valid_for<"
                                                  + detail::type_name_as_string<seq_legal_alph_type>
                                                  + "> evaluated to false on " + detail::make_printable(c)};
                            }
                            return c;
                        })
                    | views::char_to<std::ranges::range_value_t<seq_type>>, // convert to actual target alphabet
                std::back_inserter(sequence));
        }
        else
        {
            detail::consume(stream_view | detail::take_until_or_throw(is_end)); // consume until "//"
        }

        std::ranges::advance(stream_it, 3u, std::ranges::end(stream_view)); // Skip `//` and potentially '\n'
    }

    //!\copydoc sequence_file_output_format::write_sequence_record
    template <typename stream_type, // constraints checked by file
              typename seq_type,    // other constraints checked inside function
              typename id_type,
              typename qual_type>
    void write_sequence_record(stream_type & stream,
                               sequence_file_output_options const & options,
                               seq_type && sequence,
                               id_type && id,
                               qual_type && SEQAN3_DOXYGEN_ONLY(qualities))
    {
        std::ostreambuf_iterator stream_it{stream};
        size_t sequence_size{0};
        [[maybe_unused]] char buffer[50];
        if constexpr (!detail::decays_to_ignore_v<seq_type>)
            sequence_size = std::ranges::size(sequence);

        // ID
        if constexpr (detail::decays_to_ignore_v<id_type>)
        {
            throw std::logic_error{"The ID field may not be set to ignore when writing genbank files."};
        }
        else if (std::ranges::empty(id)) //[[unlikely]]
        {
            throw std::runtime_error{"The ID field may not be empty when writing genbank files."};
        }
        else if (options.embl_genbank_complete_header)
        {
            std::ranges::copy(id, stream_it);
        }
        else
        {
            std::ranges::copy(std::string_view{"LOCUS       "}, stream_it);
            std::ranges::copy(id, stream_it);
            std::ranges::copy(std::string_view{"                 "}, stream_it);
            auto res = std::to_chars(&buffer[0], &buffer[0] + sizeof(buffer), sequence_size);
            std::copy(&buffer[0], res.ptr, stream_it);
            std::ranges::copy(std::string_view{" bp\n"}, stream_it);
        }

        // Sequence
        if constexpr (detail::decays_to_ignore_v<seq_type>) // sequence
        {
            throw std::logic_error{"The SEQ field may not be set to ignore when writing genbank files."};
        }
        else if (std::ranges::empty(sequence)) //[[unlikely]]
        {
            throw std::runtime_error{"The SEQ field may not be empty when writing genbank files."};
        }
        else
        {
            std::ranges::copy(std::string_view{"ORIGIN\n"}, stream_it);
            auto seq = sequence | seqan3::views::chunk(60);
            size_t i = 0;
            size_t bp = 1;

            while (bp < sequence_size)
            {
                // Sequence length with more than 9 digits are not possible in one genbank entry, maximal 350 kb are
                // allowed. See: https://www.ncbi.nlm.nih.gov/Sitemap/samplerecord.html#SequenceLengthA
                for (size_t j = std::to_string(bp).size(); j < 9; j++)
                    stream_it = ' ';
                std::ranges::copy(std::to_string(bp), stream_it);
                stream_it = ' ';
                std::ranges::copy(seq[i] | views::to_char | views::interleave(10, std::string_view{" "}), stream_it);
                bp += 60;
                ++i;
                detail::write_eol(stream_it, false);
            }
            std::ranges::copy(std::string_view{"//"}, stream_it);
            detail::write_eol(stream_it, false);
        }
    }
};

} // namespace seqan3
