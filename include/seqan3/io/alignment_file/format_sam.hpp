// -----------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2019, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2019, Knut Reinert & MPI für molekulare Genetik
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at: https://github.com/seqan/seqan3/blob/master/LICENSE
// -----------------------------------------------------------------------------------------------------

/*!\file
 * \brief Provides the seqan3::alignment_file_format_sam class.
 * \author Svenja Mehringer <avenja.mehringer AT fu-berlin.de>
 */

#pragma once

#include <iterator>
#include <string>
#include <string_view>
#include <vector>

#include <range/v3/algorithm/copy.hpp>
#include <range/v3/view/remove_if.hpp>

#include <seqan3/core/concept/tuple.hpp>
#include <seqan3/core/metafunction/range.hpp>
#include <seqan3/io/alignment_file/detail.hpp>
#include <seqan3/io/alignment_file/header.hpp>
// #include <seqan3/io/alignment/input_options.hpp>
#include <seqan3/io/alignment_file/output_options.hpp>
#include <seqan3/io/alignment_file/sam_tag_dictionary.hpp>
#include <seqan3/io/detail/ignore_output_iterator.hpp>
#include <seqan3/io/detail/misc.hpp>
#include <seqan3/io/stream/parse_condition.hpp>
#include <seqan3/range/view/char_to.hpp>
#include <seqan3/range/view/get.hpp>
#include <seqan3/range/view/take_until.hpp>
#include <seqan3/range/view/to_char.hpp>
#include <seqan3/std/ranges>
#include <seqan3/std/concepts>

namespace seqan3
{

/*!\brief       The SAM format.
 * \implements  alignment_file_format_concept
 * \ingroup     alignment_file
 *
 * \details
 *
 * ### Introduction
 *
 * SAM is often used for storing alignments of several read sequences against one
 * or more reference sequences. See the
 * [article on wikipedia](https://en.wikipedia.org/wiki/SAM_(file_format)) for an
 * introduction of the format or look into the official
 * [SAM format specifications](https://samtools.github.io/hts-specs/SAMv1.pdf).
 * **SeqAn implements version 1.6 of the SAM specification**.
 *
 * ### Fields
 *
 * The SAM format provides the following fields:
 * seqan3::field::ALIGNMENT, seqan3::field::SEQ, seqan3::field::QUAL,
 * seqan3::field::ID, seqan3::field::REF_SEQ, seqan3::field::REF_ID
 * seqan3::field::REF_OSSFET, seqan3::field::OFFSET, seqan3::field::FLAG,
 * seqan3::field::MAPQ and seqan3::field::MATE.
 * In addition there is the seqan3::field::HEADER_PTR, which is usually not set but
 * needed to provide the range-based functionality of the file.
 *
 * **None of the fields are required** when writing but will be defaulted
 * to '0' for numeric fields and '*' for other fields.
 *
 * ### SAM format columns -> fields
 *
 * As many users will be accustomed to the columns of the SAM format, here is a
 * mapping of the common SAM format columns to the SeqAn3 record fields:
 *
 * | #  | SAM Column ID |  FIELD name                                       |
 * |:--:|:--------------|:--------------------------------------------------|
 * | 1  | QNAME         | seqan3::field::ID                                 |
 * | 2  | FLAG          | seqan3::field::FLAG                               |
 * | 3  | RNAME         | seqan3::field::REF_ID                             |
 * | 4  | POS           | seqan3::field::REF_OFFSET                         |
 * | 5  | MAPQ          | seqan3::field::MAPQ                               |
 * | 6  | CIGAR         | implicilty stored in seqan3::field::ALIGNMENT     |
 * | 7  | RNEXT         | seqan3::field::MATE (tuple pos 0)                 |
 * | 8  | PNEXT         | seqan3::field::MATE (tuple pos 1)                 |
 * | 9  | TLEN          | seqan3::field::MATE (tuple pos 2)                 |
 * | 10 | SEQ           | seqan3::field::SEQ                                |
 * | 11 | QUAL          | seqan3::field::QUAL                               |
 *
 * The (read sequence/query) **OFFSET** will be required to store the soft
 * clipping information at the read start (end clipping will be automatically
 * deduced by how much the read sequence length + offset is larger than the
 * alignment length).
 *
 * Note: SeqAn currently does not support hard clipping. When reading SAM,
 * hard-clipping is discarded; but the resulting alignment/sequence combination
 * is still valid.
 *
 * ### Format Check
 *
 * The format checks are implemented according to the official
 * [SAM format specifications](https://samtools.github.io/hts-specs/SAMv1.pdf)
 * in order to ensure correct SAM file output.
 *
 * If a non-recoverable format violation is encountered on reading, or you specify
 * invalid values/combinations when writing, seqan3::format_error is thrown.
 *
 * Note: All sequence like fields in SAM (e.g. field::SEQ) are truncated at the
 *       the first white space character (see seqan3::is_space) to ensure a
 *       correct format.
 *
 * ### Header implementation
 *
 * The SAM header is printed once in the beginning, before the first record is
 * written.
 */
class alignment_file_format_sam
{
public:
    /*!\name Constructors, destructor and assignment
     * \{
     */
    alignment_file_format_sam() = default;                                         //!< Defaulted
    //!\brief Copy construction is explicitly deleted, because you can't have multiple access to the same file.
    alignment_file_format_sam(alignment_file_format_sam const &) = delete;
    //!\brief Copy assignment is explicitly deleted, because you can't have multiple access to the same file.
    alignment_file_format_sam & operator=(alignment_file_format_sam const &) = delete;
    alignment_file_format_sam(alignment_file_format_sam &&) = default;             //!< Defaulted
    alignment_file_format_sam & operator=(alignment_file_format_sam &&) = default; //!< Defaulted
    ~alignment_file_format_sam() = default;                                        //!< Defaulted
    //!\}

    //!\brief The valid file extensions for this format; note that you can modify this value.
    static inline std::vector<std::string> file_extensions
    {
        { "sam" },
    };

    //!\copydoc AlignmentFileOutputFormat::write
    template <typename stream_type,
              typename seq_type,
              typename id_type,
              typename offset_type,
              typename ref_seq_type,
              typename ref_id_type,
              typename ref_offset_type,
              typename align_type,
              typename flag_type,
              typename mapq_type,
              typename qual_type,
              typename mate_type,
              typename tag_dict_type,
              typename e_value_type,
              typename bit_score_type>
    void write(stream_type                            &  stream,
               alignment_file_output_options const    &  options,
               std::unique_ptr<alignment_file_header> &  header_ptr,
               seq_type                               && seq,
               qual_type                              && qual,
               id_type                                && id,
               offset_type                            && offset,
               ref_seq_type                           && SEQAN3_DOXYGEN_ONLY(ref_seq),
               ref_id_type                            && ref_id,
               ref_offset_type                        && ref_offset,
               align_type                             && align,
               flag_type                              && flag,
               mapq_type                              && mapq,
               mate_type                              && mate,
               tag_dict_type                          && tag_dict,
               e_value_type                           && SEQAN3_DOXYGEN_ONLY(e_value),
               bit_score_type                         && SEQAN3_DOXYGEN_ONLY(bit_score))
    {
        /* Note the following general things:
         *
         * - Given the SAM specifications, all fields may be empty
         *
         * - Arithmetic values default to 0 while all others default to '*'
         *
         * - Because of the former, arithmetic values can be directly streamed
         *   into 'stream' as operator<< is defined for all arithmetic types
         *   and the default value (0) is also the SAM default.
         *
         * - All other non-arithmetic values need to be checked for emptiness
         */

        // ---------------------------------------------------------------------
        // Type Requirements (as static asserts for user friendliness)
        // ---------------------------------------------------------------------
        static_assert((std::ranges::ForwardRange<seq_type>        &&
                      Alphabet<value_type_t<remove_cvref_t<seq_type>>>),
                      "The seq object must be a std::ranges::ForwardRange over "
                      "letters that model seqan3::Alphabet.");

        static_assert((std::ranges::ForwardRange<id_type>         &&
                      Alphabet<value_type_t<remove_cvref_t<id_type>>>),
                      "The id object must be a std::ranges::ForwardRange over "
                      "letters that model seqan3::Alphabet.");

        static_assert(std::UnsignedIntegral<remove_cvref_t<offset_type>>,
                      "The offset object must be a std::UnsignedIntegral.");

        static_assert((std::ranges::ForwardRange<ref_seq_type>    &&
                      Alphabet<value_type_t<remove_cvref_t<ref_seq_type>>>),
                      "The ref_seq object must be a std::ranges::ForwardRange "
                      "over letters that model seqan3::Alphabet.");

        static_assert((std::ranges::ForwardRange<ref_id_type>     &&
                      Alphabet<value_type_t<remove_cvref_t<ref_id_type>>>),
                      "The ref_id object must be a std::ranges::ForwardRange "
                      "over letters that model seqan3::Alphabet.");

        static_assert(std::Integral<remove_cvref_t<ref_offset_type>>, // -1 is given default to evaluate to 0
                      "The ref_offset object must be an std::Integral >= 0.");

        if (((ref_offset + 1) < 0))
            throw format_error("The ref_offset object must be an std::Integral >= 0.");

        static_assert(tuple_like_concept<remove_cvref_t<align_type>>,
                      "The align object must be a std::pair of two ranges whose "
                      "value_type is comparable to seqan3::gap");

        static_assert((std::tuple_size_v<remove_cvref_t<align_type>> == 2 &&
                       std::EqualityComparableWith<gap, value_type_t<remove_cvref_t<decltype(std::get<0>(align))>>> &&
                       std::EqualityComparableWith<gap, value_type_t<remove_cvref_t<decltype(std::get<1>(align))>>>),
                      "The align object must be a std::pair of two ranges whose "
                      "value_type is comparable to seqan3::gap");

        static_assert(std::UnsignedIntegral<remove_cvref_t<flag_type>>,
                      "The flag object must be a std::UnsignedIntegral.");

        static_assert(std::UnsignedIntegral<remove_cvref_t<mapq_type>>,
                      "The mapq object must be a std::UnsignedIntegral.");

        static_assert((std::ranges::ForwardRange<qual_type>       &&
                       Alphabet<value_type_t<remove_cvref_t<qual_type>>>),
                      "The qual object must be a std::ranges::ForwardRange "
                      "over letters that model seqan3::Alphabet.");

        static_assert(tuple_like_concept<remove_cvref_t<mate_type>>,
                      "The mate object must be a std::tuple of size 3 with "
                      "1) a std::ranges::ForwardRange with a value_type modelling seqan3::Alphabet, "
                      "2) an std::UnsignedIntegral, and"
                      "3) an std::UnsignedIntegral.");

        static_assert((std::ranges::ForwardRange<decltype(std::get<0>(mate))>                     &&
                      Alphabet<value_type_t<remove_cvref_t<decltype(std::get<0>(mate))>>> &&
                      std::UnsignedIntegral<remove_cvref_t<decltype(std::get<1>(mate))>>          &&
                      std::UnsignedIntegral<remove_cvref_t<decltype(std::get<2>(mate))>>),
                      "The mate object must be a std::tuple of size 3 with "
                      "1) a std::ranges::ForwardRange with a value_type modelling seqan3::Alphabet, "
                      "2) an std::UnsignedIntegral, and"
                      "3) an std::UnsignedIntegral.");

        static_assert(std::Same<remove_cvref_t<tag_dict_type>, sam_tag_dictionary>,
                      "The tag_dict object must be of type seqan3::sam_tag_dictionary.");

        // ---------------------------------------------------------------------
        // logical Requirements
        // ---------------------------------------------------------------------
        if (!empty(get<1>(align)) && empty(seq))
            throw format_error("If you specify an align object you must also specify the seq object. "
                               "Hint: Check if offset needs to be set to if soft-clipping is present.");

        if (options.sam_require_header && (header_ptr != nullptr) && !empty(ref_id))
        {
            if ((header_ptr->ref_dict).count(std::string(ref_id)) == 0) // no reference id matched
                throw format_error(std::string("The ref_id '") + std::string(ref_id) +
                                   "' was not in the list of references");
        }

        // ---------------------------------------------------------------------
        // Writing the Header on first call
        // ---------------------------------------------------------------------
        if (options.sam_require_header && !written_header && (header_ptr != nullptr))
        {
            write_header(stream, options, header_ptr);
            written_header = true;
        }

        // ---------------------------------------------------------------------
        // Writing the Record
        // ---------------------------------------------------------------------
        std::ranges::ostreambuf_iterator stream_it{stream};
        char const separator{'\t'};

        write_range(stream_it, std::forward<id_type>(id));

        stream << separator;

        stream << std::forward<flag_type>(flag) << separator;

        write_range(stream_it, std::forward<ref_id_type>(ref_id));

        stream << separator;

        stream << (ref_offset + 1) << separator; // SAM is 1 based

        stream << std::forward<mapq_type>(mapq) << separator;

        if (!empty(get<1>(align)))
        {
            // compute possible distance from alignment end to sequence end
            // which indicates soft clipping at the end.
            // This should be replace by a free count_gaps function for
            // aligned sequences which is more efficient if possible.
            size_t off_end{seq.size() - offset};
            for (auto chr : get<1>(align))
                if (chr == gap{})
                    ++off_end;
            off_end -= (get<1>(align)).size();

            write_range(stream_it,
                        detail::get_cigar_string(std::forward<align_type>(align),
                                                 std::forward<offset_type>(offset),
                                                 off_end));
        }
        else
        {
            stream << '*';
        }

        stream << separator;

        write_range(stream_it, get<0>(std::forward<mate_type>(mate)));

        stream << separator;

        stream << get<1>(std::forward<mate_type>(mate)) << separator;

        stream << get<2>(std::forward<mate_type>(mate)) << separator;

        write_range(stream_it, std::forward<seq_type>(seq));

        stream << separator;

        write_range(stream_it, std::forward<qual_type>(qual));

        write_tag_fields(stream, std::forward<tag_dict_type>(tag_dict), separator);

        detail::write_eol(stream_it, options.add_carriage_return);
    }

protected:
    //!\privatesection
    //!\brief The format version string.
    static constexpr char format_version[4] = "1.6";

    //!\brief A variable that tracks whether the header_ptr as been written or not.
    bool written_header{false};

    /*!\brief Writes a field value to the stream.
     * \tparam stream_it_t The stream iterator type.
     * \tparam field_type  The type of the field value. Must model std::ranges::ForwardRange.
     *
     * \param[in,out] stream_it   The stream iterator to print to.
     * \param[in]     field_value The value to print.
     */
    template <typename stream_it_t, typename field_type>
    //!\cond
        requires std::ranges::ForwardRange<field_type>
    //!\endcond
    void write_range(stream_it_t & stream_it, field_type && field_value)
    {
        if (empty(field_value))
            stream_it = '*';
        else
            std::ranges::copy(field_value | view::to_char | view::take_until(is_space), stream_it);
    }

    /*!\brief Writes the optional fields of the seqan3::sam_tag_dictionary.
     * \tparam stream_t   The stream type.
     *
     * \param[in,out] stream    The stream to print to.
     * \param[in]     tag_dict  The tag dictionary to print.
     * \param[in]     separator The field separator to append.
     */
    template <typename stream_t>
    void write_tag_fields(stream_t & stream, sam_tag_dictionary const & tag_dict, char const separator)
    {
        auto stream_variant_fn = [&stream] (auto && arg) // helper to print an std::variant
        {
            using T = remove_cvref_t<decltype(arg)>;

            if constexpr (!container_concept<T>)
            {
                stream << arg;
            }
            else
            {
                if (arg.begin() != arg.end())
                {
                    for (auto it = arg.begin(); it != (arg.end() - 1); ++it)
                        stream << *it << ",";

                    stream << *(arg.end() - 1); // write last value without trailing ','
                }
            }
        };

        for (auto & [tag, variant] : tag_dict)
        {
            stream << separator;

            char char0 = tag / 256;
            char char1 = tag % 256;

            stream << char0 << char1 << ':' << detail::sam_tag_type_char[variant.index()] << ':';

            if (detail::sam_tag_type_char_extra[variant.index()] != '\0')
                stream << detail::sam_tag_type_char_extra[variant.index()] << ',';

            std::visit(stream_variant_fn, variant);
        }
    }

    /*!\brief Writes the SAM header_ptr.
     * \tparam stream_t   The stream type.
     *
     * \param[in,out] stream  The stream to print to.
     * \param[in]     options The options to alter printing.
     * \param[in]     header_ptr  The header_ptr (as a pointer) to print.
     *
     * \throws seqan3::format_error if the header_ptr object contains the wrong
     *         information or the contents are ill-formed.
     *
     * \details
     *
     * Before writing the header_ptr, the contents are checked for correctness
     * according to the rules of the official
     * [SAM format specifications](https://samtools.github.io/hts-specs/SAMv1.pdf).
     */
    template <typename stream_t>
    void write_header(stream_t & stream,
                      alignment_file_output_options const & options,
                      std::unique_ptr<alignment_file_header> & header_ptr)
    {
        if (header_ptr != nullptr)
        {
            // -----------------------------------------------------------------
            // Check Header
            // -----------------------------------------------------------------

            // (@HD) Check header_ptr line
            // The format version string will be taken from the local member variable
            if (!(header_ptr->sorting == "unknown"   ||
                  header_ptr->sorting == "unsorted"  ||
                  header_ptr->sorting == "queryname" ||
                  header_ptr->sorting == "coordinate" ))
                throw format_error{"SAM format error: The header_ptr->sorting member must be "
                                   "one of [unknown, unsorted, queryname, coordinate]."};

            if (!(header_ptr->grouping == "none"   ||
                  header_ptr->grouping == "query"  ||
                  header_ptr->grouping == "reference"))
                throw format_error{"SAM format error: The header_ptr->grouping member must be "
                                   "one of [none, query, reference]."};

            // (@SQ) Check Reference Sequence Dictionary lines

            // TODO

            // - sorting order be one of ...
            // - grouping can be one of ...
            // - reference names must be unique
            // - ids of read groups must be unique
            // - program ids need to be unique
            // many more small semantic things, like fits REGEX

            // -----------------------------------------------------------------
            // Write Header
            // -----------------------------------------------------------------
            std::ranges::ostreambuf_iterator stream_it{stream};

            // (@HD) Write header_ptr line [required].
            stream << "@HD\tVN:";
            stream << format_version;

            if (!header_ptr->sorting.empty())
                stream << "\tSO:" << header_ptr->sorting;

            if (!header_ptr->grouping.empty())
                stream << "\tGO:" << header_ptr->grouping;

            detail::write_eol(stream_it, options.add_carriage_return);

            // (@SQ) Write Reference Sequence Dictionary lines [required].
            for (auto const & [ref_name, ref_info] : header_ptr->ref_dict)
            {
                stream << "@SQ"
                       << "\tSN:" << ref_name
                       << "\tLN:" << get<0>(ref_info);

                if (!get<1>(ref_info).empty())
                    stream << "\t" << get<1>(ref_info);

                detail::write_eol(stream_it, options.add_carriage_return);
            }

            // Write read group (@RG) lines if specified.
            for (auto const & read_group : header_ptr->read_groups)
            {
                stream << "@RG"
                       << "\tID:" << get<0>(read_group);

                if (!get<1>(read_group).empty())
                    stream << "\t" << get<1>(read_group);

                detail::write_eol(stream_it, options.add_carriage_return);
            }

            // Write program (@PG) lines if specified.
            for (auto const & program : header_ptr->program_infos)
            {
                stream << "@PG"
                       << "\tID:" << program.id;

                if (!program.name.empty())
                    stream << "\tPN:" << program.name;

                if (!program.command_line_call.empty())
                    stream << "\tCL:" << program.command_line_call;

                if (!program.previous.empty())
                    stream << "\tPP:" << program.previous;

                if (!program.description.empty())
                    stream << "\tDS:" << program.description;

                if (!program.version.empty())
                    stream << "\tVN:" << program.version;

                detail::write_eol(stream_it, options.add_carriage_return);
            }

            // Write comment (@CO) lines if specified.
            for (auto const & comment : header_ptr->comments)
            {
                stream << "@CO\t" << comment;
                detail::write_eol(stream_it, options.add_carriage_return);
            }
        }
    }
};

} // namespace seqan3
