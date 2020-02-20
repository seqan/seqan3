// -----------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2020, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2020, Knut Reinert & MPI für molekulare Genetik
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at: https://github.com/seqan/seqan3/blob/master/LICENSE.md
// -----------------------------------------------------------------------------------------------------

/*!\file
 * \brief Provides the seqan3::sequence_file_format_genbank class.
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
#include <seqan3/range/views/interleave.hpp>
#include <seqan3/range/views/istreambuf.hpp>
#include <seqan3/range/views/to_char.hpp>
#include <seqan3/range/views/take.hpp>
#include <seqan3/range/views/take_line.hpp>
#include <seqan3/range/views/take_until.hpp>
#include <seqan3/std/algorithm>
#include <seqan3/std/charconv>
#include <seqan3/std/ranges>

namespace seqan3
{

/*!\brief       The GenBank format.
 * \implements  seqan3::sequence_file_input_format
 * \implements  seqan3::sequence_file_output_format
 * \ingroup     sequence
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
 */
class format_genbank
{
public:
    /*!\name Constructors, destructor and assignment
     * \{
     */
    format_genbank() noexcept = default; //!< Defaulted.
    format_genbank(format_genbank const &) noexcept = default; //!< Defaulted.
    format_genbank & operator=(format_genbank const &) noexcept = default; //!< Defaulted.
    format_genbank(format_genbank &&) noexcept = default; //!< Defaulted.
    format_genbank & operator=(format_genbank &&) noexcept = default; //!< Defaulted.
    ~format_genbank() noexcept = default; //!< Defaulted.
    //!\}

    //!\brief The valid file extensions for this format; note that you can modify this value.
    static inline std::vector<std::string> file_extensions
    {
        { "genbank" },
        { "gb" },
        { "gbk" },
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

        if (!(std::ranges::equal(stream_view | views::take_until_or_throw(is_cntrl || is_blank), std::string{"LOCUS"})))
            throw parse_error{"An entry has to start with the code word LOCUS."};

        //ID
        if constexpr (!detail::decays_to_ignore_v<id_type>)
        {
            if (options.embl_genbank_complete_header)
            {
                std::ranges::copy(std::string_view{"LOCUS"}, std::ranges::back_inserter(id));

                while (!is_char<'O'>(*std::ranges::begin(stream_view)))
                {
                        std::ranges::copy(stream_view | views::take_line_or_throw
                                                      | views::char_to<std::ranges::range_value_t<id_type>>,
                                                        std::ranges::back_inserter(id));
                        id.push_back('\n');
                }
            }
            else
            {
                detail::consume(stream_view | views::take_until(!is_blank));

                auto read_id_until = [&stream_view, &id] (auto predicate)
                {
                    std::ranges::copy(stream_view | views::take_until_or_throw(predicate)
                                                  | views::char_to<std::ranges::range_value_t<id_type>>,
                                      std::ranges::back_inserter(id));
                };

                if (options.truncate_ids)
                    read_id_until(is_space);
                else
                    read_id_until(is_cntrl);

                detail::consume(stream_view | views::take_line_or_throw);
            }
        }

        // Jump to sequence
        while (!(is_char<'O'>(*std::ranges::begin(stream_view)) || options.embl_genbank_complete_header))
            detail::consume(stream_view | views::take_line_or_throw);

        // Sequence
        detail::consume(stream_view | views::take_line_or_throw); // consume "ORIGIN"
        auto constexpr is_end = is_char<'/'> ;
        if constexpr (!detail::decays_to_ignore_v<seq_type>)
        {
            auto constexpr is_legal_alph = is_in_alphabet<seq_legal_alph_type>;
            std::ranges::copy(stream_view | std::views::filter(!(is_space || is_digit))
                                          | views::take_until_or_throw_and_consume(is_end) // consume "//"
                                          | std::views::transform([is_legal_alph] (char const c) // enforce legal alphabet
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
                                          | views::char_to<std::ranges::range_value_t<seq_type>>,    // convert to actual target alphabet
                                            std::ranges::back_inserter(sequence));
        }
        else
        {
            detail::consume(stream_view | views::take_until_or_throw_and_consume(is_end)); // consume until "//"
            ++stream_it; // consume "/n"
        }
    }

    //!\copydoc sequence_file_output_format::write_sequence_record
    template <typename stream_type,     // constraints checked by file
              typename seq_type,        // other constraints checked inside function
              typename id_type,
              typename qual_type>
    void write_sequence_record(stream_type                        & stream,
                               sequence_file_output_options const & options,
                               seq_type                           && sequence,
                               id_type                            && id,
                               qual_type                          && SEQAN3_DOXYGEN_ONLY(qualities))
    {
        seqan3::ostreambuf_iterator stream_it{stream};
        size_t sequence_size{0};
        [[maybe_unused]] char buffer[50];
        if constexpr (!detail::decays_to_ignore_v<seq_type>)
            sequence_size = ranges::size(sequence);

        // ID
        if constexpr (detail::decays_to_ignore_v<id_type>)
        {
            throw std::logic_error{"The ID field may not be set to ignore when writing genbank files."};
        }
        else if (ranges::empty(id)) //[[unlikely]]
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
            auto seq = sequence | ranges::view::chunk(60);
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
                std::ranges::copy(seq[i] | views::to_char
                                         | views::interleave(10, std::string_view{" "}), stream_it);
                bp += 60;
                ++i;
                detail::write_eol(stream_it,false);
            }
            std::ranges::copy(std::string_view{"//"}, stream_it);
            detail::write_eol(stream_it,false);
        }
    }
};

} // namespace seqan
