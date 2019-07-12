// -----------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2019, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2019, Knut Reinert & MPI für molekulare Genetik
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at: https://github.com/seqan/seqan3/blob/master/LICENSE.md
// -----------------------------------------------------------------------------------------------------

/*!\file
 * \brief Provides the seqan3::format_vienna tag and the seqan3::sequencestrucure_file_input_format and
 *        seqan3::strucure_file_output_format specialisation for this tag.
 * \author Jörg Winkler <j.winkler AT fu-berlin.de>
 */

#pragma once

#include <cstdio>
#include <iterator>
#include <stack>
#include <string>
#include <string_view>
#include <type_traits>
#include <vector>

#include <range/v3/algorithm/copy.hpp>
#include <range/v3/view/chunk.hpp>
#include <range/v3/view/drop_while.hpp>
#include <range/v3/view/join.hpp>
#include <range/v3/view/remove_if.hpp>
#include <range/v3/view/transform.hpp>

#include <seqan3/alphabet/structure/wuss.hpp>
#include <seqan3/core/char_operations/predicate.hpp>
#include <seqan3/core/type_traits/range.hpp>
#include <seqan3/io/detail/misc.hpp>
#include <seqan3/io/stream/iterator.hpp>
#include <seqan3/io/structure_file/detail.hpp>
#include <seqan3/io/structure_file/input_format_concept.hpp>
#include <seqan3/io/structure_file/input_options.hpp>
#include <seqan3/io/structure_file/output_format_concept.hpp>
#include <seqan3/io/structure_file/output_options.hpp>
#include <seqan3/range/shortcuts.hpp>
#include <seqan3/range/detail/misc.hpp>
#include <seqan3/range/view/char_to.hpp>
#include <seqan3/range/view/istreambuf.hpp>
#include <seqan3/range/view/to_char.hpp>
#include <seqan3/range/view/take.hpp>
#include <seqan3/range/view/take_line.hpp>
#include <seqan3/range/view/take_until.hpp>
#include <seqan3/std/algorithm>
#include <seqan3/std/ranges>

namespace seqan3
{
/*!\brief       The Vienna format (dot bracket notation) for RNA sequences with secondary structure.
 * \implements  seqan3::StructureFileInputFormat
 * \implements  seqan3::StructureFileOutputFormat
 * \ingroup     structure_file
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
 * The Vienna format provides the fields seqan3::field::SEQ, seqan3::field::ID, seqan3::field::BPP (read only),
 * seqan3::field::STRUCTURE, seqan3::field::STRUCTURED_SEQ and seqan3::field::ENERGY.\n\n
 * If you select seqan3::field::STRUCTURED_SEQ you must not select seqan3::field::SEQ or seqan3::field::STRUCTURE.\n
 * Either the field seqan3::field::SEQ or the field seqan3::field::STRUCTURED_SEQ is required when writing.\n
 * The field seqan3::field::BPP is ignored when writing, but a structure string can be converted to bpp when reading.
 *
 * ### Implementation notes
 *
 * When reading the ID-line the identifier (`>`) and any blank characters before the actual ID are
 * stripped. Each field is read/written as a single line (except ENERGY, which goes right after the structure).
 * Numbers and spaces within the sequence are simply ignored, but not within the structure.
 */
struct format_vienna
{
    //!\brief The valid file extensions for this format; note that you can modify this value.
    static inline std::vector<std::string> file_extensions
    {
        { "dbn" },
        { "fasta" },
        { "fa" }
    };
};

} // namespace seqan

namespace seqan3::detail
{

//!\brief The seqan3::structure_file_input_format specialisation that can handle formatted VIENNA input.
//!\ingroup structure_file
template <>
class structure_file_input_format<format_vienna>
{
public:
    //!\brief Exposes the format tag that this class is specialised with.
    using format_tag = format_vienna;

    /*!\name Constructors, destructor and assignment
     * \{
     */
    structure_file_input_format()                                                noexcept = default; //!< Defaulted.
    //!\brief Copy construction is explicitly deleted, because you can't have multiple access to the same file.
    structure_file_input_format(structure_file_input_format const &)                      = delete;
    //!\brief Copy assignment is explicitly deleted, because you can't have multiple access to the same file.
    structure_file_input_format & operator=(structure_file_input_format const &)          = delete;
    structure_file_input_format(structure_file_input_format &&)                  noexcept = default; //!< Defaulted.
    structure_file_input_format & operator=(structure_file_input_format &&)      noexcept = default; //!< Defaulted.
    ~structure_file_input_format()                                               noexcept = default; //!< Defaulted.
    //!\}

    //!\copydoc seqan3::StructureFileInputFormat::read
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
        auto stream_view = view::istreambuf(stream);

        // READ ID (if present)
        auto constexpr is_id = is_char<'>'>;
        if (is_id(*begin(stream_view)))
        {
            if constexpr (!detail::decays_to_ignore_v<id_type>)
            {
                if (options.truncate_ids)
                {
                    std::ranges::copy(stream_view | std::view::drop_while(is_id || is_blank) // skip leading >
                                                  | view::take_until_or_throw(is_cntrl || is_blank)
                                                  | view::char_to<value_type_t<id_type>>,
                                      std::back_inserter(id));
                    detail::consume(stream_view | view::take_line_or_throw);
                }
                else
                {
                    std::ranges::copy(stream_view | std::view::drop_while(is_id || is_blank) // skip leading >
                                                  | view::take_line_or_throw
                                                  | view::char_to<value_type_t<id_type>>,
                                      std::back_inserter(id));
                }
            }
            else
            {
                detail::consume(stream_view | view::take_line_or_throw);
            }
        }
        else if constexpr (!detail::decays_to_ignore_v<id_type>)
        {
            auto constexpr is_legal_seq = is_in_alphabet<seq_legal_alph_type>;
            if (!is_legal_seq(*begin(stream_view))) // if neither id nor seq found: throw
            {
                throw parse_error{std::string{"Expected to be on beginning of ID or sequence, but "} +
                                  is_id.msg.str() + " and " + is_legal_seq.msg.str() +
                                  " evaluated to false on " + detail::make_printable(*begin(stream_view))};
            }
        }

        // READ SEQUENCE
        if constexpr (!detail::decays_to_ignore_v<seq_type>)
        {
            auto constexpr is_legal_seq = is_in_alphabet<seq_legal_alph_type>;
            std::ranges::copy(stream_view | view::take_line_or_throw                      // until end of line
                                          | std::view::filter(!(is_space || is_digit)) // ignore whitespace and numbers
                                          | std::view::transform([is_legal_seq](char const c)
                                            {
                                                if (!is_legal_seq(c))                    // enforce legal alphabet
                                                {
                                                    throw parse_error{std::string{"Encountered an unexpected letter: "} +
                                                                      is_legal_seq.msg.str() +
                                                                      " evaluated to false on " +
                                                                      detail::make_printable(c)};
                                                }
                                              return c;
                                            })
                                          | view::char_to<value_type_t<seq_type>>, // convert to actual target alphabet
                              std::back_inserter(seq));
        }
        else
        {
            detail::consume(stream_view | view::take_line_or_throw);
        }

        // READ STRUCTURE (if present)
        if constexpr (!detail::decays_to_ignore_v<structure_type>)
        {
            if constexpr (structured_seq_combined)
            {
                assert(std::addressof(seq) == std::addressof(structure));
                using alph_type = typename value_type_t<structure_type>::structure_alphabet_type;
                std::ranges::copy(read_structure<alph_type>(stream_view), begin(structure));

                if constexpr (!detail::decays_to_ignore_v<bpp_type>)
                    detail::bpp_from_rna_structure<alph_type>(bpp, structure);
            }
            else
            {
                using alph_type = value_type_t<structure_type>;
                std::ranges::copy(read_structure<alph_type>(stream_view), std::back_inserter(structure));

                if constexpr (!detail::decays_to_ignore_v<bpp_type>)
                    detail::bpp_from_rna_structure<alph_type>(bpp, structure);
            }
            if constexpr (!detail::decays_to_ignore_v<seq_type>)
                if (size(seq) != size(structure))
                    throw parse_error{"Found sequence and associated structure of different length."};
        }
        else if constexpr (!detail::decays_to_ignore_v<bpp_type>)
        {
            detail::bpp_from_rna_structure<wuss51>(bpp, read_structure<wuss51>(stream_view));

            if constexpr (!detail::decays_to_ignore_v<seq_type>)
                if (size(seq) != size(bpp))
                    throw parse_error{"Found sequence and associated structure of different length."};
        }
        else
        {
            detail::consume(stream_view | view::take_until(is_space)); // until whitespace
        }

        // READ ENERGY (if present)
        if constexpr (!detail::decays_to_ignore_v<energy_type>)
        {
            std::string e_str = stream_view | view::take_line
                                            | std::view::filter(!(is_space || is_char<'('> || is_char<')'>));
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
        }
        detail::consume(stream_view | view::take_until(!is_space));
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
        auto constexpr is_legal_structure = is_in_alphabet<alph_type>;
        return stream_view | view::take_until(is_space) // until whitespace
                           | std::view::transform([is_legal_structure](char const c)
                             {
                                 if (!is_legal_structure(c))
                                 {
                                     throw parse_error{
                                         std::string{"Encountered an unexpected letter: "} +
                                         is_legal_structure.msg.str() +
                                         " evaluated to false on " + detail::make_printable(c)};
                                 }
                                 return c;
                             })                                  // enforce legal alphabet
                           | view::char_to<alph_type>;           // convert to actual target alphabet
    }
};

//!\brief The seqan3::structure_file_output_format specialisation that can write formatted VIENNA.
//!\ingroup structure_file
template <>
class structure_file_output_format<format_vienna>
{
public:
    //!\brief Exposes the format tag that this class is specialised with.
    using format_tag = format_vienna;

    /*!\name Constructors, destructor and assignment
     * \{
     */
    structure_file_output_format()                                                 noexcept = default; //!< Defaulted.
    //!\brief Copy construction is explicitly deleted, because you can't have multiple access to the same file.
    structure_file_output_format(structure_file_output_format const &)                      = delete;
    //!\brief Copy assignment is explicitly deleted, because you can't have multiple access to the same file.
    structure_file_output_format & operator=(structure_file_output_format const &)          = delete;
    structure_file_output_format(structure_file_output_format &&)                  noexcept = default; //!< Defaulted.
    structure_file_output_format & operator=(structure_file_output_format &&)      noexcept = default; //!< Defaulted.
    ~structure_file_output_format()                                                noexcept = default; //!< Defaulted.
    //!\}

    //!\copydoc seqan3::StructureFileOutputFormat::write
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
        seqan3::ostreambuf_iterator stream_it{stream};

        // WRITE ID (optional)
        if constexpr (!detail::decays_to_ignore_v<id_type>)
        {
            if (!empty(id))
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
            if (empty(seq)) //[[unlikely]]
                throw std::runtime_error{"The SEQ field may not be empty when writing Vienna files."};

            std::ranges::copy(seq | view::to_char, stream_it);
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
            if (!empty(structure))
                std::ranges::copy(structure | view::to_char, stream_it);

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
};

} // namespace seqan3::detail
