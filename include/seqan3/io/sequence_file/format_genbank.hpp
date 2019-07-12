// -----------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2019, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2019, Knut Reinert & MPI für molekulare Genetik
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
#include <seqan3/range/view/char_to.hpp>
#include <seqan3/range/view/interleave.hpp>
#include <seqan3/range/view/istreambuf.hpp>
#include <seqan3/range/view/to_char.hpp>
#include <seqan3/range/view/take.hpp>
#include <seqan3/range/view/take_line.hpp>
#include <seqan3/range/view/take_until.hpp>
#include <seqan3/std/algorithm>
#include <seqan3/std/charconv>
#include <seqan3/std/ranges>

namespace seqan3
{

/*!\brief       The GenBank format (tag).
 * \implements  seqan3::SequenceFileInputFormat
 * \implements  seqan3::SequenceFileOutputFormat
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
 * There is no truncate_ids option while reading because the GenBank format has no (FASTA-like) idbuffer. Instead,
 * there is the option "complete_header" to indicate whether the whole header is to be read
 * (embl_genbank_complete_header=true) or only the "LOCUS" information should be stored.
 *
 * Qualities passed to the write function are ignored.
 *
 */
struct format_genbank
{
    //!\brief The valid file extensions for this format; note that you can modify this value.
    static inline std::vector<std::string> file_extensions
    {
        { "genbank" },
        { "gb" },
        { "gbk" },
    };
};

} // namespace seqan

namespace seqan3::detail
{

//!\brief The seqan3::sequence_file_input_format specialisation that handles formatted Genbank input.
//!\ingroup sequence
template<>
class sequence_file_input_format<format_genbank>
{
public:
    //!\brief Exposes the format tag that this class is specialised with.
    using format_tag = format_genbank;

    /*!\name Constructors, destructor and assignment
     * \{
     */
    sequence_file_input_format()                                               noexcept = default; //!< Defaulted.
    //!\brief Copy construction is explicitly deleted, because you can't have multiple access to the same file.
    sequence_file_input_format(sequence_file_input_format const &)                      = delete;
    //!\brief Copy assignment is explicitly deleted, because you can't have multiple access to the same file.
    sequence_file_input_format & operator=(sequence_file_input_format const &)          = delete;
    sequence_file_input_format(sequence_file_input_format &&)                  noexcept = default; //!< Defaulted.
    sequence_file_input_format & operator=(sequence_file_input_format &&)      noexcept = default; //!< Defaulted.
    ~sequence_file_input_format()                                              noexcept = default; //!< Defaulted.
    //!\}

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
              qual_type                                                                 & SEQAN3_DOXYGEN_ONLY(qualities))
    {
        auto stream_view = view::istreambuf(stream);
        auto stream_it = std::ranges::begin(stream_view);

        if (!(std::ranges::equal(stream_view | view::take_until_or_throw(is_cntrl || is_blank), std::string{"LOCUS"})))
            throw parse_error{"An entry has to start with the code word LOCUS."};

        //ID
        if constexpr (!detail::decays_to_ignore_v<id_type>)
        {
            if (options.embl_genbank_complete_header)
            {
                std::ranges::copy(std::string_view{"LOCUS"}, std::back_inserter(id));

                while (!is_char<'O'>(*std::ranges::begin(stream_view)))
                {
                        std::ranges::copy(stream_view | view::take_line_or_throw
                                                      | view::char_to<value_type_t<id_type>>,
                                                        std::back_inserter(id));
                        id.push_back('\n');
                }
            }
            else
            {
                detail::consume(stream_view | view::take_until(!is_blank));
                // read id
                std::ranges::copy(stream_view | view::take_until_or_throw(is_cntrl || is_blank)
                                              | view::char_to<value_type_t<id_type>>,
                                                std::back_inserter(id));
                detail::consume(stream_view | view::take_line_or_throw);
            }
        }

        // Jump to sequence
        while (!(is_char<'O'>(*std::ranges::begin(stream_view)) || options.embl_genbank_complete_header))
            detail::consume(stream_view | view::take_line_or_throw);

        // Sequence
        detail::consume(stream_view | view::take_line_or_throw); // consume "ORIGIN"
        auto constexpr is_end = is_char<'/'> ;
        if constexpr (!detail::decays_to_ignore_v<seq_type>)
        {
            auto constexpr is_legal_alph = is_in_alphabet<seq_legal_alph_type>;
            std::ranges::copy(stream_view | std::view::filter(!(is_space || is_digit))
                                          | view::take_until_or_throw_and_consume(is_end) // consume "//"
                                          | std::view::transform([is_legal_alph] (char const c) // enforce legal alphabet
                                            {
                                                if (!is_legal_alph(c))
                                                {
                                                    throw parse_error{std::string{"Encountered an unexpected letter: "} +
                                                                      is_legal_alph.msg.str() +
                                                                      " evaluated to false on " +
                                                                      detail::make_printable(c)};
                                                }
                                                return c;
                                            })
                                          | view::char_to<value_type_t<seq_type>>,    // convert to actual target alphabet
                                            std::back_inserter(sequence));
        }
        else
        {
            detail::consume(stream_view | view::take_until_or_throw_and_consume(is_end)); // consume until "//"
            ++stream_it; // consume "/n"
        }
    }
};

//!\brief The seqan3::sequence_file_output_format specialisation that can write formatted Genbank output.
//!\ingroup sequence
template <>
class sequence_file_output_format<format_genbank>
{
public:
    //!\brief Exposes the format tag that this class is specialised with.
    using format_tag = format_genbank;

    /*!\name Constructors, destructor and assignment
     * \{
     */
    sequence_file_output_format()                                                noexcept = default; //!< Defaulted.
    //!\brief Copy construction is explicitly deleted, because you can't have multiple access to the same file.
    sequence_file_output_format(sequence_file_output_format const &)                      = delete;
    //!\brief Copy assignment is explicitly deleted, because you can't have multiple access to the same file.
    sequence_file_output_format & operator=(sequence_file_output_format const &)          = delete;
    sequence_file_output_format(sequence_file_output_format &&)                  noexcept = default; //!< Defaulted.
    sequence_file_output_format & operator=(sequence_file_output_format &&)      noexcept = default; //!< Defaulted.
    ~sequence_file_output_format()                                               noexcept = default; //!< Defaulted.
    //!\}

    //!\copydoc SequenceFileOutputFormat::write
    template <typename stream_type,     // constraints checked by file
              typename seq_type,        // other constraints checked inside function
              typename id_type,
              typename qual_type>
    void write(stream_type                        & stream,
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
                std::ranges::copy(seq[i] | view::to_char
                                         | view::interleave(10, std::string_view{" "}), stream_it);
                bp += 60;
                ++i;
                detail::write_eol(stream_it,false);
            }
            std::ranges::copy(std::string_view{"//"}, stream_it);
            detail::write_eol(stream_it,false);
        }
    }
};

} // namespace seqan3::detail
