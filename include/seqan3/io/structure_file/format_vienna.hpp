// SPDX-FileCopyrightText: 2006-2025 Knut Reinert & Freie Universität Berlin
// SPDX-FileCopyrightText: 2016-2025 Knut Reinert & MPI für molekulare Genetik
// SPDX-License-Identifier: BSD-3-Clause

/*!\file
 * \brief Provides the seqan3::format_vienna.
 * \author Jörg Winkler <j.winkler AT fu-berlin.de>
 */

#pragma once

#include <algorithm>
#include <cstdio>
#include <iterator>
#include <ranges>
#include <stack>
#include <string>
#include <string_view>
#include <type_traits>
#include <vector>

#include <seqan3/alphabet/adaptation/char.hpp>
#include <seqan3/alphabet/concept.hpp>
#include <seqan3/alphabet/structure/wuss.hpp>
#include <seqan3/alphabet/views/char_to.hpp>
#include <seqan3/alphabet/views/to_char.hpp>
#include <seqan3/core/range/detail/misc.hpp>
#include <seqan3/core/range/type_traits.hpp>
#include <seqan3/io/detail/misc.hpp>
#include <seqan3/io/structure_file/detail.hpp>
#include <seqan3/io/structure_file/input_format_concept.hpp>
#include <seqan3/io/structure_file/input_options.hpp>
#include <seqan3/io/structure_file/output_format_concept.hpp>
#include <seqan3/io/structure_file/output_options.hpp>
#include <seqan3/io/views/detail/istreambuf_view.hpp>
#include <seqan3/io/views/detail/take_line_view.hpp>
#include <seqan3/io/views/detail/take_until_view.hpp>
#include <seqan3/utility/char_operations/predicate.hpp>
#include <seqan3/utility/detail/type_name_as_string.hpp>
#include <seqan3/utility/range/to.hpp>

namespace seqan3
{
/*!\brief       The Vienna format (dot bracket notation) for RNA sequences with secondary structure.
 * \implements  seqan3::structure_file_input_format
 * \implements  seqan3::structure_file_output_format
 * \ingroup     io_structure_file
 *
 * \details
 *
 * ### Introduction
 *
 * Dot Bracket or Vienna Notation is widely used for secondary structure annotation. Is is a very simple format,
 * containing one or more sequences. Each sequence must appear as a single line in the file.
 * A sequence may be preceded by a special line starting with the '>' character followed by a sequence name
 * (like FASTA). After each sequence line there is usually a line containing
 * secondary structure, using brackets to denote interacting nucleotides or amino acids, and dots for unpaired sites.
 * The length of the struture must equal the length of the sequence.
 * Optionally, the structure may be followed by a space character and the minimum free energy value enclosed
 * in parentheses (). Note that there cannot be energy without structure.
 *
 * The Vienna format is the output format of _RNAfold_. Furthermore, it is designed to be compatible with the
 * input format of the ViennaRNA package (if structure and energy are omitted).
 * See https://www.tbi.univie.ac.at/RNA/tutorial/#sec2_7 for details.
 *
 * ### fields_specialisation
 *
 * The Vienna format provides the fields seqan3::field::seq, seqan3::field::id, seqan3::field::bpp (read only),
 * seqan3::field::structure, seqan3::field::structured_seq and seqan3::field::energy.\n\n
 * If you select seqan3::field::structured_seq you must not select seqan3::field::seq or seqan3::field::structure.\n
 * Either the field seqan3::field::seq or the field seqan3::field::structured_seq is required when writing.\n
 * The field seqan3::field::bpp is ignored when writing, but a structure string can be converted to bpp when reading.
 *
 * ### Implementation notes
 *
 * When reading the ID-line the identifier (`>`) and any blank characters before the actual ID are
 * stripped. Each field is read/written as a single line (except ENERGY, which goes right after the structure).
 * Numbers and spaces within the sequence are simply ignored, but not within the structure.
 *
 * \remark For a complete overview, take a look at \ref io_structure_file
 */
class format_vienna
{
public:
    /*!\name Constructors, destructor and assignment
     * \{
     */
    format_vienna() noexcept = default;                                  //!< Defaulted.
    format_vienna(format_vienna const &) noexcept = default;             //!< Defaulted.
    format_vienna & operator=(format_vienna const &) noexcept = default; //!< Defaulted.
    format_vienna(format_vienna &&) noexcept = default;                  //!< Defaulted.
    format_vienna & operator=(format_vienna &&) noexcept = default;      //!< Defaulted.
    ~format_vienna() noexcept = default;                                 //!< Defaulted.

    //!\}

    //!\brief The valid file extensions for this format; note that you can modify this value.
    static inline std::vector<std::string> file_extensions{{"dbn"}, {"fasta"}, {"fa"}};

protected:
    //!\copydoc seqan3::structure_file_input_format::read_structure_record
    template <typename stream_type, // constraints checked by file
              typename seq_legal_alph_type,
              bool structured_seq_combined,
              typename seq_type, // other constraints checked inside function
              typename id_type,
              typename bpp_type,
              typename structure_type,
              typename energy_type,
              typename react_type,
              typename comment_type,
              typename offset_type>
    void
    read_structure_record(stream_type & stream,
                          structure_file_input_options<seq_legal_alph_type, structured_seq_combined> const & options,
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
        auto stream_view = detail::istreambuf(stream);

        // READ ID (if present)
        constexpr auto is_id = is_char<'>'>;
        if (is_id(*begin(stream_view)))
        {
            if constexpr (!detail::decays_to_ignore_v<id_type>)
            {
                if (options.truncate_ids)
                {
                    std::ranges::copy(stream_view | std::views::drop_while(is_id || is_blank) // skip leading >
                                          | detail::take_until_or_throw(is_cntrl || is_blank)
                                          | views::char_to<std::ranges::range_value_t<id_type>>,
                                      std::back_inserter(id));
                    detail::consume(stream_view | detail::take_line_or_throw);
                }
                else
                {
                    std::ranges::copy(stream_view | std::views::drop_while(is_id || is_blank) // skip leading >
                                          | detail::take_line_or_throw
                                          | views::char_to<std::ranges::range_value_t<id_type>>,
                                      std::back_inserter(id));
                }
            }
            else
            {
                detail::consume(stream_view | detail::take_line_or_throw);
            }
        }
        else if constexpr (!detail::decays_to_ignore_v<id_type>)
        {
            constexpr auto is_legal_seq = char_is_valid_for<seq_legal_alph_type>;
            if (!is_legal_seq(*begin(stream_view))) // if neither id nor seq found: throw
            {
                throw parse_error{std::string{"Expected to be on beginning of ID or sequence, but "} + is_id.msg
                                  + " and char_is_valid_for<" + detail::type_name_as_string<seq_legal_alph_type>
                                  + ">" + " evaluated to false on " + detail::make_printable(*begin(stream_view))};
            }
        }

        // READ SEQUENCE
        if constexpr (!detail::decays_to_ignore_v<seq_type>)
        {
            constexpr auto is_legal_seq = char_is_valid_for<seq_legal_alph_type>;
            std::ranges::copy(
                stream_view | detail::take_line_or_throw          // until end of line
                    | std::views::filter(!(is_space || is_digit)) // ignore whitespace and numbers
                    | std::views::transform(
                        [is_legal_seq](char const c)
                        {
                            if (!is_legal_seq(c)) // enforce legal alphabet
                            {
                                throw parse_error{std::string{"Encountered an unexpected letter: "}
                                                  + "char_is_valid_for<"
                                                  + detail::type_name_as_string<seq_legal_alph_type>
                                                  + "> evaluated to false on " + detail::make_printable(c)};
                            }
                            return c;
                        })
                    | views::char_to<std::ranges::range_value_t<seq_type>>, // convert to actual target alphabet
                std::back_inserter(seq));
        }
        else
        {
            detail::consume(stream_view | detail::take_line_or_throw);
        }

        // READ STRUCTURE (if present)
        [[maybe_unused]] int64_t structure_length{};
        if constexpr (!detail::decays_to_ignore_v<structure_type>)
        {
            if constexpr (structured_seq_combined)
            {
                assert(std::addressof(seq) == std::addressof(structure));
                using alph_type = typename std::ranges::range_value_t<structure_type>::structure_alphabet_type;
                // We need the structure_length parameter to count the length of the structure while reading
                // because we cannot infer it from the (already resized) structure_seq object.
                auto range = read_structure<alph_type>(stream_view);
                // Use std::views::take to avoid going out of bounds if the structure is longer than the sequence.
                auto res = std::ranges::copy(range | std::views::take(std::ranges::distance(seq)),
                                             std::ranges::begin(structure));
                structure_length = std::ranges::distance(std::ranges::begin(structure), res.out);
                // If the structure is longer than the sequence, there are characters left.
                // std::ranges::distance will also consume the characters in the stream.
                structure_length += std::ranges::distance(range);

                if constexpr (!detail::decays_to_ignore_v<bpp_type>)
                    detail::bpp_from_rna_structure<alph_type>(bpp, structure);
            }
            else
            {
                using alph_type = std::ranges::range_value_t<structure_type>;
                std::ranges::copy(read_structure<alph_type>(stream_view), std::back_inserter(structure));
                structure_length = std::ranges::distance(structure);

                if constexpr (!detail::decays_to_ignore_v<bpp_type>)
                    detail::bpp_from_rna_structure<alph_type>(bpp, structure);
            }
        }
        else if constexpr (!detail::decays_to_ignore_v<bpp_type>)
        {
            detail::bpp_from_rna_structure<wuss51>(bpp, read_structure<wuss51>(stream_view));
            structure_length = std::ranges::distance(bpp);
        }
        else
        {
            detail::consume(stream_view | detail::take_until(is_space)); // until whitespace
        }

        if constexpr (!detail::decays_to_ignore_v<seq_type>
                      && !(detail::decays_to_ignore_v<structure_type> && detail::decays_to_ignore_v<bpp_type>))
        {
            if (std::ranges::distance(seq) != structure_length)
                throw parse_error{"Found sequence and associated structure of different length."};
        }

        // READ ENERGY (if present)
        if constexpr (!detail::decays_to_ignore_v<energy_type>)
        {
            std::string e_str = stream_view | detail::take_line
                              | std::views::filter(!(is_space || is_char<'('> || is_char<')'>))
                              | ranges::to<std::string>();
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
            detail::consume(stream_view | detail::take_line);
        }
        detail::consume(stream_view | detail::take_until(!is_space));
    }

    //!\copydoc seqan3::structure_file_output_format::write_structure_record
    template <typename stream_type, // constraints checked by file
              typename seq_type,    // other constraints checked inside function
              typename id_type,
              typename bpp_type,
              typename structure_type,
              typename energy_type,
              typename react_type,
              typename comment_type,
              typename offset_type>
    void write_structure_record(stream_type & stream,
                                structure_file_output_options const & options,
                                seq_type && seq,
                                id_type && id,
                                bpp_type && SEQAN3_DOXYGEN_ONLY(bpp),
                                structure_type && structure,
                                energy_type && energy,
                                react_type && SEQAN3_DOXYGEN_ONLY(react),
                                react_type && SEQAN3_DOXYGEN_ONLY(react_err),
                                comment_type && SEQAN3_DOXYGEN_ONLY(comment),
                                offset_type && SEQAN3_DOXYGEN_ONLY(offset))
    {
        std::ostreambuf_iterator stream_it{stream};

        // WRITE ID (optional)
        if constexpr (!detail::decays_to_ignore_v<id_type>)
        {
            if (!std::ranges::empty(id))
            {
                stream_it = '>';
                stream_it = ' ';
                std::ranges::copy(id, stream_it);
                detail::write_eol(stream_it, options.add_carriage_return);
            }
        }

        // WRITE SEQUENCE
        if constexpr (!detail::decays_to_ignore_v<seq_type>)
        {
            if (std::ranges::empty(seq)) //[[unlikely]]
                throw std::runtime_error{"The SEQ field may not be empty when writing Vienna files."};

            std::ranges::copy(seq | views::to_char, stream_it);
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
            if (!std::ranges::empty(structure))
                std::ranges::copy(structure | views::to_char, stream_it);

            // WRITE ENERGY (optional)
            if constexpr (!detail::decays_to_ignore_v<energy_type>)
            {
                if (energy)
                {
                    // TODO(joergi-w) enable the following when std::to_chars is implemented for float types
                    //                    auto [endptr, ec] = std::to_chars(str.data(),
                    //                                                      str.data() + str.size(),
                    //                                                      energy,
                    //                                                      std::chars_format::fixed,
                    //                                                      options.precision);
                    //                    if (ec == std::errc())
                    //                        std::ranges::copy(str.data(), endptr, stream_it);
                    //                    else
                    //                        throw std::runtime_error{"The energy could not be transformed into a string."};

                    stream_it = ' ';
                    stream_it = '(';

                    std::array<char, 100> str;
                    int len = std::snprintf(str.data(), 100, "%.*f", options.precision, energy);
                    if (len < 0 || len >= 100)
                        throw std::runtime_error{"The energy could not be transformed into a string."};
                    std::ranges::copy(str.data(), str.data() + len, stream_it);

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
    /*!\brief Extract the structure string from the given stream.
     * \tparam alph_type        The alphabet type the structure is converted to.
     * \tparam stream_view_type The type of the input stream.
     * \param stream_view       The input stream to be read.
     * \return                  A std::ranges::view containing the structure annotation string.
     */
    template <typename alph_type, typename stream_view_type>
    auto read_structure(stream_view_type & stream_view)
    {
        constexpr auto is_legal_structure = char_is_valid_for<alph_type>;
        return stream_view | detail::take_until(is_space) // until whitespace
             | std::views::transform(
                   [is_legal_structure](char const c)
                   {
                       if (!is_legal_structure(c))
                       {
                           throw parse_error{std::string{"Encountered an unexpected letter: char_is_valid_for<"}
                                             + detail::type_name_as_string<alph_type>
                                             + "> evaluated to false on " + detail::make_printable(c)};
                       }
                       return c;
                   })                     // enforce legal alphabet
             | views::char_to<alph_type>; // convert to actual target alphabet
    }
};

} // namespace seqan3
