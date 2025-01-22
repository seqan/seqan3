// SPDX-FileCopyrightText: 2006-2025 Knut Reinert & Freie Universität Berlin
// SPDX-FileCopyrightText: 2016-2025 Knut Reinert & MPI für molekulare Genetik
// SPDX-License-Identifier: BSD-3-Clause

/*!\file
 * \brief Provides the seqan3::format_sam.
 * \author Svenja Mehringer <svenja.mehringer AT fu-berlin.de>
 */

#pragma once

#include <iterator>
#include <ranges>
#include <string>
#include <vector>

#include <seqan3/alphabet/concept.hpp>
#include <seqan3/alphabet/views/to_char.hpp>
#include <seqan3/io/sam_file/detail/cigar.hpp>
#include <seqan3/io/sam_file/detail/format_sam_base.hpp>
#include <seqan3/io/sam_file/header.hpp>
#include <seqan3/io/sam_file/input_format_concept.hpp>
#include <seqan3/io/sam_file/output_options.hpp>
#include <seqan3/io/sam_file/sam_flag.hpp>
#include <seqan3/io/sam_file/sam_tag_dictionary.hpp>
#include <seqan3/io/sequence_file/input_format_concept.hpp>
#include <seqan3/io/sequence_file/output_options.hpp>
#include <seqan3/io/stream/detail/fast_ostreambuf_iterator.hpp>
#include <seqan3/io/views/detail/istreambuf_view.hpp>
#include <seqan3/io/views/detail/take_until_view.hpp>
#include <seqan3/utility/detail/type_name_as_string.hpp>
#include <seqan3/utility/range/to.hpp>
#include <seqan3/utility/tuple/concept.hpp>
#include <seqan3/utility/views/slice.hpp>

namespace seqan3
{

/*!\brief       The SAM format (tag).
 * \implements  AlignmentFileFormat
 * \ingroup io_sam_file
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
 * Take a look at our tutorial \ref tutorial_sam_file for a walk through of how to read SAM/BAM files.
 *
 * ### fields_specialisation
 *
 * The SAM format provides the following fields:
 * seqan3::field::seq, seqan3::field::qual,
 * seqan3::field::id, seqan3::field::ref_seq, seqan3::field::ref_id
 * seqan3::field::ref_offset, seqan3::field::offset, seqan3::field::flag,
 * seqan3::field::mapq and seqan3::field::mate.
 * In addition there is the seqan3::field::header_ptr, which is usually only used internally
 * to provide the range-based functionality of the file.
 *
 * **None of the fields are required** when writing. If they are not given, a default value of '0' for numeric fields
 * and '*' for other fields is used.
 *
 * ### SAM format columns -> fields
 *
 * Since many users will be accustomed to the columns of the SAM format, here is a
 * mapping of the common SAM format columns to the SeqAn record fields:
 *
 * | #  | SAM Column ID |  FIELD name                       |
 * |:--:|:--------------|:----------------------------------|
 * | 1  | QNAME         | seqan3::field::id                 |
 * | 2  | FLAG          | seqan3::field::flag               |
 * | 3  | RNAME         | seqan3::field::ref_id             |
 * | 4  | POS           | seqan3::field::ref_offset         |
 * | 5  | MAPQ          | seqan3::field::mapq               |
 * | 6  | CIGAR         | seqan3::field::cigar              |
 * | 7  | RNEXT         | seqan3::field::mate (tuple pos 0) |
 * | 8  | PNEXT         | seqan3::field::mate (tuple pos 1) |
 * | 9  | TLEN          | seqan3::field::mate (tuple pos 2) |
 * | 10 | SEQ           | seqan3::field::seq                |
 * | 11 | QUAL          | seqan3::field::qual               |
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
 * ### Header implementation
 *
 * The SAM header (if present) is read/written once in the beginning before the
 * first record is read/written.
 *
 * \remark For a complete overview, take a look at \ref io_sam_file
 */
class format_sam : protected detail::format_sam_base
{
public:
    /*!\name Constructors, destructor and assignment
     * \{
     */
    // construction cannot be noexcept because this class has a std::string variable as a quality string buffer.
    format_sam() = default;                              //!< Defaulted.
    format_sam(format_sam const &) = delete;             //!< Deleted. Header holds a unique_ptr.
    format_sam & operator=(format_sam const &) = delete; //!< Deleted. Header holds a unique_ptr.
    format_sam(format_sam &&) = default;                 //!< Defaulted.
    format_sam & operator=(format_sam &&) = default;     //!< Defaulted.
    ~format_sam() = default;                             //!< Defaulted.

    //!\}

    //!\brief The valid file extensions for this format; note that you can modify this value.
    static inline std::vector<std::string> file_extensions{
        {"sam"},
    };

protected:
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
                              qual_type & qualities);

    template <typename stream_type, // constraints checked by file
              typename seq_type,    // other constraints checked inside function
              typename id_type,
              typename qual_type>
    void write_sequence_record(stream_type & stream,
                               sequence_file_output_options const & SEQAN3_DOXYGEN_ONLY(options),
                               seq_type && sequence,
                               id_type && id,
                               qual_type && qualities);

    template <typename stream_type, // constraints checked by file
              typename seq_legal_alph_type,
              typename ref_seqs_type,
              typename ref_ids_type,
              typename stream_pos_type,
              typename seq_type,
              typename id_type,
              typename ref_seq_type,
              typename ref_id_type,
              typename ref_offset_type,
              typename cigar_type,
              typename flag_type,
              typename mapq_type,
              typename qual_type,
              typename mate_type,
              typename tag_dict_type,
              typename e_value_type,
              typename bit_score_type>
    void read_alignment_record(stream_type & stream,
                               sam_file_input_options<seq_legal_alph_type> const & SEQAN3_DOXYGEN_ONLY(options),
                               ref_seqs_type & ref_seqs,
                               sam_file_header<ref_ids_type> & header,
                               stream_pos_type & position_buffer,
                               seq_type & seq,
                               qual_type & qual,
                               id_type & id,
                               ref_seq_type & SEQAN3_DOXYGEN_ONLY(ref_seq),
                               ref_id_type & ref_id,
                               ref_offset_type & ref_offset,
                               cigar_type & cigar_vector,
                               flag_type & flag,
                               mapq_type & mapq,
                               mate_type & mate,
                               tag_dict_type & tag_dict,
                               e_value_type & SEQAN3_DOXYGEN_ONLY(e_value),
                               bit_score_type & SEQAN3_DOXYGEN_ONLY(bit_score));

    template <typename stream_type,
              typename header_type,
              typename seq_type,
              typename id_type,
              typename ref_seq_type,
              typename ref_id_type,
              typename qual_type,
              typename mate_type,
              typename tag_dict_type,
              typename e_value_type,
              typename bit_score_type>
    void write_alignment_record(stream_type & stream,
                                sam_file_output_options const & options,
                                header_type && header,
                                seq_type && seq,
                                qual_type && qual,
                                id_type && id,
                                ref_seq_type && SEQAN3_DOXYGEN_ONLY(ref_seq),
                                ref_id_type && ref_id,
                                std::optional<int32_t> ref_offset,
                                std::vector<cigar> const & cigar_vector,
                                sam_flag const flag,
                                uint8_t const mapq,
                                mate_type && mate,
                                tag_dict_type && tag_dict,
                                e_value_type && SEQAN3_DOXYGEN_ONLY(e_value),
                                bit_score_type && SEQAN3_DOXYGEN_ONLY(bit_score));

private:
    //!\brief Stores quality values temporarily if seq and qual information are combined (not supported by SAM yet).
    std::string tmp_qual{};

    //!\brief An empty dummy container to pass to align_format.write() such that an empty field is written.
    static constexpr std::string_view dummy{};

    //!\brief The default header for the alignment format.
    sam_file_header<> default_header{};

    //!\brief A buffer to store a raw record pointing into the stream buffer of the input.
    std::array<std::string_view, 11> raw_record{};

    //!brief Returns a reference to dummy if passed a std::ignore.
    std::string_view const & default_or(detail::ignore_t) const noexcept
    {
        return dummy;
    }

    //!brief Returns the input unchanged.
    template <typename t>
    decltype(auto) default_or(t && v) const noexcept
    {
        return std::forward<t>(v);
    }

    template <arithmetic value_type>
    void read_sam_dict_vector(seqan3::detail::sam_tag_variant & variant, std::string_view const str, value_type value);

    void read_sam_byte_vector(seqan3::detail::sam_tag_variant & variant, std::string_view const str);

    void read_sam_dict(std::string_view const tag_str, sam_tag_dictionary & target);

    template <typename stream_it_t, std::ranges::forward_range field_type>
    void write_range_or_asterisk(stream_it_t & stream_it, field_type && field_value);

    template <typename stream_it_t>
    void write_range_or_asterisk(stream_it_t & stream_it, char const * const field_value);

    template <typename stream_it_t>
    void write_tag_fields(stream_it_t & stream, sam_tag_dictionary const & tag_dict, char const separator);
};

//!\copydoc sequence_file_input_format::read_sequence_record
template <typename stream_type, // constraints checked by file
          typename seq_legal_alph_type,
          typename stream_pos_type,
          typename seq_type, // other constraints checked inside function
          typename id_type,
          typename qual_type>
inline void format_sam::read_sequence_record(stream_type & stream,
                                             sequence_file_input_options<seq_legal_alph_type> const & options,
                                             stream_pos_type & position_buffer,
                                             seq_type & sequence,
                                             id_type & id,
                                             qual_type & qualities)
{
    sam_file_input_options<seq_legal_alph_type> align_options;

    {
        read_alignment_record(stream,
                              align_options,
                              std::ignore,
                              default_header,
                              position_buffer,
                              sequence,
                              qualities,
                              id,
                              std::ignore,
                              std::ignore,
                              std::ignore,
                              std::ignore,
                              std::ignore,
                              std::ignore,
                              std::ignore,
                              std::ignore,
                              std::ignore,
                              std::ignore);
    }

    if constexpr (!detail::decays_to_ignore_v<seq_type>)
        if (std::ranges::distance(sequence) == 0)
            throw parse_error{"The sequence information must not be empty."};
    if constexpr (!detail::decays_to_ignore_v<id_type>)
    {
        if (std::ranges::distance(id) == 0)
            throw parse_error{"The id information must not be empty."};
        if (options.truncate_ids)
            id = id | detail::take_until_and_consume(is_space) | ranges::to<id_type>();
    }
}

//!\copydoc sequence_file_output_format::write_sequence_record
template <typename stream_type, // constraints checked by file
          typename seq_type,    // other constraints checked inside function
          typename id_type,
          typename qual_type>
inline void format_sam::write_sequence_record(stream_type & stream,
                                              sequence_file_output_options const & SEQAN3_DOXYGEN_ONLY(options),
                                              seq_type && sequence,
                                              id_type && id,
                                              qual_type && qualities)
{
    using default_mate_t = std::tuple<std::string_view, std::optional<int32_t>, int32_t>;

    sam_file_output_options output_options;

    write_alignment_record(stream,
                           output_options,
                           /*header*/ std::ignore,
                           /*seq*/ default_or(sequence),
                           /*qual*/ default_or(qualities),
                           /*id*/ default_or(id),
                           /*ref_seq*/ std::string_view{},
                           /*ref_id*/ std::string_view{},
                           /*ref_offset*/ -1,
                           /*cigar_vector*/ std::vector<cigar>{},
                           /*flag*/ sam_flag::none,
                           /*mapq*/ 0,
                           /*mate*/ default_mate_t{},
                           /*tag_dict*/ sam_tag_dictionary{},
                           /*e_value*/ 0,
                           /*bit_score*/ 0);
}

//!\copydoc sam_file_input_format::read_alignment_record
template <typename stream_type, // constraints checked by file
          typename seq_legal_alph_type,
          typename ref_seqs_type,
          typename ref_ids_type,
          typename stream_pos_type,
          typename seq_type,
          typename id_type,
          typename ref_seq_type,
          typename ref_id_type,
          typename ref_offset_type,
          typename cigar_type,
          typename flag_type,
          typename mapq_type,
          typename qual_type,
          typename mate_type,
          typename tag_dict_type,
          typename e_value_type,
          typename bit_score_type>
inline void
format_sam::read_alignment_record(stream_type & stream,
                                  sam_file_input_options<seq_legal_alph_type> const & SEQAN3_DOXYGEN_ONLY(options),
                                  ref_seqs_type & ref_seqs,
                                  sam_file_header<ref_ids_type> & header,
                                  stream_pos_type & position_buffer,
                                  seq_type & seq,
                                  qual_type & qual,
                                  id_type & id,
                                  ref_seq_type & SEQAN3_DOXYGEN_ONLY(ref_seq),
                                  ref_id_type & ref_id,
                                  ref_offset_type & ref_offset,
                                  cigar_type & cigar_vector,
                                  flag_type & flag,
                                  mapq_type & mapq,
                                  mate_type & mate,
                                  tag_dict_type & tag_dict,
                                  e_value_type & SEQAN3_DOXYGEN_ONLY(e_value),
                                  bit_score_type & SEQAN3_DOXYGEN_ONLY(bit_score))
{
    static_assert(detail::decays_to_ignore_v<ref_offset_type>
                      || detail::is_type_specialisation_of_v<ref_offset_type, std::optional>,
                  "The ref_offset must be a specialisation of std::optional.");

    auto stream_it = detail::fast_istreambuf_iterator{*stream.rdbuf()};

    auto stream_view = detail::istreambuf(stream);

    int32_t ref_offset_tmp{}; // needed to read the ref_offset (int) beofre storing it in std::optional<uint32_t>
    std::ranges::range_value_t<decltype(header.ref_ids())> ref_id_tmp{}; // in case mate is requested but ref_offset not

    // Header
    // -------------------------------------------------------------------------------------------------------------
    if (is_char<'@'>(*stream_it)) // we always read the header if present
    {
        read_header(stream_view, header, ref_seqs);

        if (std::ranges::begin(stream_view) == std::ranges::end(stream_view)) // file has no records
            return;
    }

    // Store the current file position in the buffer.
    position_buffer = stream.tellg();

    // We don't know wether we have 11 or 12 fields, since the tags are optional.
    // The last field will thus contain either the quality sequence
    // or the quality sequence AND tags. This will be handled at the respective fields below.
    stream_it.cache_record_into('\n', '\t', raw_record);

    // Fields 1-5: ID FLAG REF_ID REF_OFFSET MAPQ
    // -------------------------------------------------------------------------------------------------------------
    if constexpr (!detail::decays_to_ignore_v<id_type>)
        read_forward_range_field(raw_record[0], id);

    uint16_t flag_integral{};
    read_arithmetic_field(raw_record[1], flag_integral);
    flag = sam_flag{flag_integral};

    read_forward_range_field(raw_record[2], ref_id_tmp);
    check_and_assign_ref_id(ref_id, ref_id_tmp, header, ref_seqs);

    read_arithmetic_field(raw_record[3], ref_offset_tmp);
    --ref_offset_tmp; // SAM format is 1-based but SeqAn operates 0-based

    if (ref_offset_tmp == -1)
        ref_offset = std::nullopt; // indicates an unmapped read -> ref_offset is not set
    else if (ref_offset_tmp > -1)
        ref_offset = ref_offset_tmp;
    else if (ref_offset_tmp < -1)
        throw format_error{"No negative values are allowed for field::ref_offset."};

    if constexpr (!detail::decays_to_ignore_v<mapq_type>)
        read_arithmetic_field(raw_record[4], mapq);

    // Field 6: CIGAR
    // -------------------------------------------------------------------------------------------------------------
    if constexpr (!detail::decays_to_ignore_v<cigar_type>)
        cigar_vector = detail::parse_cigar(raw_record[5]);

    // Field 7-9: (RNEXT PNEXT TLEN) = MATE
    // -------------------------------------------------------------------------------------------------------------
    if constexpr (!detail::decays_to_ignore_v<mate_type>)
    {
        std::ranges::range_value_t<decltype(header.ref_ids())> tmp_mate_ref_id{};
        read_forward_range_field(raw_record[6], tmp_mate_ref_id); // RNEXT

        if (tmp_mate_ref_id == "=") // indicates "same as ref id"
        {
            if constexpr (!detail::decays_to_ignore_v<ref_id_type>)
                get<0>(mate) = ref_id;
            else
                check_and_assign_ref_id(get<0>(mate), ref_id_tmp, header, ref_seqs);
        }
        else
        {
            check_and_assign_ref_id(get<0>(mate), tmp_mate_ref_id, header, ref_seqs);
        }

        int32_t tmp_pnext{};
        read_arithmetic_field(raw_record[7], tmp_pnext); // PNEXT

        if (tmp_pnext > 0)
            get<1>(mate) = --tmp_pnext; // SAM format is 1-based but SeqAn operates 0-based.
        else if (tmp_pnext < 0)
            throw format_error{"No negative values are allowed at the mate mapping position."};
        // tmp_pnext == 0 indicates an unmapped mate -> do not fill std::optional get<1>(mate)

        read_arithmetic_field(raw_record[8], get<2>(mate)); // TLEN
    }

    // Field 10: Sequence
    // -------------------------------------------------------------------------------------------------------------
    if constexpr (!detail::decays_to_ignore_v<seq_type>)
    {
        std::string_view const seq_str = raw_record[9];

        if (!seq_str.starts_with('*')) // * indicates missing sequence information
        {
            seq.resize(seq_str.size());
            constexpr auto is_legal_alph = char_is_valid_for<seq_legal_alph_type>;

            for (size_t i = 0; i < seq_str.size(); ++i)
            {
                if (!is_legal_alph(seq_str[i]))
                    throw parse_error{std::string{"Encountered an unexpected letter: "} + "char_is_valid_for<"
                                      + detail::type_name_as_string<seq_legal_alph_type>
                                      + "> evaluated to false on " + detail::make_printable(seq_str[i])};

                seq[i] = assign_char_to(seq_str[i], std::ranges::range_value_t<seq_type>{});
            }
        }
    }

    // Field 11:  Quality
    // -------------------------------------------------------------------------------------------------------------
    // We don't know wether we have 11 or 12 fields, since the tags are optional.
    // The last field will thus contain either the quality sequence
    // or the quality sequence AND tags.
    size_t tag_begin_pos = raw_record[10].find('\t');

    std::string_view qualities =
        (tag_begin_pos == std::string_view::npos) ? raw_record[10] : raw_record[10].substr(0, tag_begin_pos);

    if constexpr (!detail::decays_to_ignore_v<qual_type>)
        read_forward_range_field(qualities, qual);

    if constexpr (!detail::decays_to_ignore_v<seq_type> && !detail::decays_to_ignore_v<qual_type>)
    {
        if (std::ranges::distance(seq) != 0 && std::ranges::distance(qual) != 0
            && std::ranges::distance(seq) != std::ranges::distance(qual))
        {
            throw format_error{detail::to_string("Sequence length (",
                                                 std::ranges::distance(seq),
                                                 ") and quality length (",
                                                 std::ranges::distance(qual),
                                                 ") must be the same.")};
        }
    }

    // All remaining optional fields if any: SAM tags dictionary
    // -------------------------------------------------------------------------------------------------------------
    if constexpr (!detail::decays_to_ignore_v<tag_dict_type>)
    {
        while (tag_begin_pos != std::string_view::npos) // read all tags if present
        {
            ++tag_begin_pos; // skip '\t'
            size_t const tag_end_pos = raw_record[10].find('\t', tag_begin_pos);

            char const * tag_begin = raw_record[10].begin() + tag_begin_pos;
            char const * tag_end =
                (tag_end_pos == std::string_view::npos) ? raw_record[10].end() : raw_record[10].begin() + tag_end_pos;

            read_sam_dict(std::string_view{tag_begin, tag_end}, tag_dict);

            tag_begin_pos = tag_end_pos;
        }
    }

    assert(stream_it == std::default_sentinel_t{} || *stream_it == '\n');
    ++stream_it; // Move from end of record to the beginning of the next or to the end of the stream.
}

//!\copydoc sam_file_output_format::write_alignment_record
template <typename stream_type,
          typename header_type,
          typename seq_type,
          typename id_type,
          typename ref_seq_type,
          typename ref_id_type,
          typename qual_type,
          typename mate_type,
          typename tag_dict_type,
          typename e_value_type,
          typename bit_score_type>
inline void format_sam::write_alignment_record(stream_type & stream,
                                               sam_file_output_options const & options,
                                               header_type && header,
                                               seq_type && seq,
                                               qual_type && qual,
                                               id_type && id,
                                               ref_seq_type && SEQAN3_DOXYGEN_ONLY(ref_seq),
                                               ref_id_type && ref_id,
                                               std::optional<int32_t> ref_offset,
                                               std::vector<cigar> const & cigar_vector,
                                               sam_flag const flag,
                                               uint8_t const mapq,
                                               mate_type && mate,
                                               tag_dict_type && tag_dict,
                                               e_value_type && SEQAN3_DOXYGEN_ONLY(e_value),
                                               bit_score_type && SEQAN3_DOXYGEN_ONLY(bit_score))
{
    /* Note the following general things:
     *
     * - Given the SAM specifications, all fields may be empty
     *
     * - arithmetic values default to 0 while all others default to '*'
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
    static_assert((std::ranges::forward_range<seq_type> && alphabet<std::ranges::range_reference_t<seq_type>>),
                  "The seq object must be a std::ranges::forward_range over "
                  "letters that model seqan3::alphabet.");

    static_assert((std::ranges::forward_range<id_type> && alphabet<std::ranges::range_reference_t<id_type>>),
                  "The id object must be a std::ranges::forward_range over "
                  "letters that model seqan3::alphabet.");

    if constexpr (!detail::decays_to_ignore_v<ref_id_type>)
    {
        static_assert((std::ranges::forward_range<ref_id_type> || std::integral<std::remove_reference_t<ref_id_type>>
                       || detail::is_type_specialisation_of_v<std::remove_cvref_t<ref_id_type>, std::optional>),
                      "The ref_id object must be a std::ranges::forward_range "
                      "over letters that model seqan3::alphabet.");

        if constexpr (std::integral<std::remove_cvref_t<ref_id_type>>
                      || detail::is_type_specialisation_of_v<std::remove_cvref_t<ref_id_type>, std::optional>)
            static_assert(!detail::decays_to_ignore_v<header_type>,
                          "If you give indices as reference id information the header must also be present.");
    }

    static_assert((std::ranges::forward_range<qual_type> && alphabet<std::ranges::range_reference_t<qual_type>>),
                  "The qual object must be a std::ranges::forward_range "
                  "over letters that model seqan3::alphabet.");

    static_assert(tuple_like<std::remove_cvref_t<mate_type>>,
                  "The mate object must be a std::tuple of size 3 with "
                  "1) a std::ranges::forward_range with a value_type modelling seqan3::alphabet, "
                  "2) a std::integral or std::optional<std::integral>, and "
                  "3) a std::integral.");

    static_assert(
        ((std::ranges::forward_range<decltype(std::get<0>(mate))>
          || std::integral<std::remove_cvref_t<decltype(std::get<0>(mate))>>
          || detail::is_type_specialisation_of_v<std::remove_cvref_t<decltype(std::get<0>(mate))>, std::optional>)
         && (std::integral<std::remove_cvref_t<decltype(std::get<1>(mate))>>
             || detail::is_type_specialisation_of_v<std::remove_cvref_t<decltype(std::get<1>(mate))>, std::optional>)
         && std::integral<std::remove_cvref_t<decltype(std::get<2>(mate))>>),
        "The mate object must be a std::tuple of size 3 with "
        "1) a std::ranges::forward_range with a value_type modelling seqan3::alphabet, "
        "2) a std::integral or std::optional<std::integral>, and "
        "3) a std::integral.");

    if constexpr (std::integral<std::remove_cvref_t<decltype(std::get<0>(mate))>>
                  || detail::is_type_specialisation_of_v<std::remove_cvref_t<decltype(std::get<0>(mate))>,
                                                         std::optional>)
        static_assert(!detail::decays_to_ignore_v<header_type>,
                      "If you give indices as mate reference id information the header must also be present.");

    static_assert(std::same_as<std::remove_cvref_t<tag_dict_type>, sam_tag_dictionary>,
                  "The tag_dict object must be of type seqan3::sam_tag_dictionary.");

    // ---------------------------------------------------------------------
    // logical Requirements
    // ---------------------------------------------------------------------
    if constexpr (!detail::decays_to_ignore_v<header_type> && !detail::decays_to_ignore_v<ref_id_type>
                  && !std::integral<std::remove_reference_t<ref_id_type>>
                  && !detail::is_type_specialisation_of_v<std::remove_reference_t<ref_id_type>, std::optional>)
    {

        if (options.sam_require_header && !std::ranges::empty(ref_id))
        {
            auto id_it = header.ref_dict.end();

            if constexpr (std::ranges::contiguous_range<decltype(ref_id)> && std::ranges::sized_range<decltype(ref_id)>
                          && std::ranges::borrowed_range<decltype(ref_id)>)
            {
                id_it = header.ref_dict.find(std::span{std::ranges::data(ref_id), std::ranges::size(ref_id)});
            }
            else
            {
                using header_ref_id_type = std::remove_reference_t<decltype(header.ref_ids()[0])>;

                static_assert(implicitly_convertible_to<ref_id_type, header_ref_id_type>,
                              "The ref_id type is not convertible to the reference id information stored in the "
                              "reference dictionary of the header object.");

                id_it = header.ref_dict.find(ref_id);
            }

            if (id_it == header.ref_dict.end()) // no reference id matched
                throw format_error{detail::to_string("The ref_id '",
                                                     ref_id,
                                                     "' was not in the list of references:",
                                                     header.ref_ids())};
        }
    }

    if (ref_offset.has_value() && (ref_offset.value() + 1) < 0)
        throw format_error{"The ref_offset object must be a std::integral >= 0."};

    // ---------------------------------------------------------------------
    // Writing the Header on first call
    // ---------------------------------------------------------------------
    if constexpr (!detail::decays_to_ignore_v<header_type>)
    {
        if (options.sam_require_header && !header_was_written)
        {
            write_header(stream, options, header);
            header_was_written = true;
        }
    }

    // ---------------------------------------------------------------------
    // Writing the Record
    // ---------------------------------------------------------------------

    detail::fast_ostreambuf_iterator stream_it{*stream.rdbuf()};
    constexpr char separator{'\t'};

    write_range_or_asterisk(stream_it, id);
    *stream_it = separator;

    stream_it.write_number(static_cast<uint16_t>(flag));
    *stream_it = separator;

    if constexpr (!detail::decays_to_ignore_v<ref_id_type>)
    {
        if constexpr (std::integral<std::remove_reference_t<ref_id_type>>)
        {
            write_range_or_asterisk(stream_it, (header.ref_ids())[ref_id]);
        }
        else if constexpr (detail::is_type_specialisation_of_v<std::remove_reference_t<ref_id_type>, std::optional>)
        {
            if (ref_id.has_value())
                write_range_or_asterisk(stream_it, (header.ref_ids())[ref_id.value()]);
            else
                *stream_it = '*';
        }
        else
        {
            write_range_or_asterisk(stream_it, ref_id);
        }
    }
    else
    {
        *stream_it = '*';
    }

    *stream_it = separator;

    // SAM is 1 based, 0 indicates unmapped read if optional is not set
    stream_it.write_number(ref_offset.value_or(-1) + 1);
    *stream_it = separator;

    stream_it.write_number(static_cast<unsigned>(mapq));
    *stream_it = separator;

    if (!std::ranges::empty(cigar_vector))
    {
        for (auto & c : cigar_vector) //TODO THIS IS PROBABLY TERRIBLE PERFORMANCE_WISE
            stream_it.write_range(c.to_string());
    }
    else
    {
        *stream_it = '*';
    }

    *stream_it = separator;

    if constexpr (std::integral<std::remove_reference_t<decltype(get<0>(mate))>>)
    {
        write_range_or_asterisk(stream_it, (header.ref_ids())[get<0>(mate)]);
    }
    else if constexpr (detail::is_type_specialisation_of_v<std::remove_reference_t<decltype(get<0>(mate))>,
                                                           std::optional>)
    {
        if (get<0>(mate).has_value())
            write_range_or_asterisk(stream_it, header.ref_ids()[get<0>(mate).value()]);
        else
            *stream_it = '*';
    }
    else
    {
        write_range_or_asterisk(stream_it, get<0>(mate));
    }

    *stream_it = separator;

    if constexpr (detail::is_type_specialisation_of_v<std::remove_cvref_t<decltype(get<1>(mate))>, std::optional>)
    {
        // SAM is 1 based, 0 indicates unmapped read if optional is not set
        stream_it.write_number(get<1>(mate).value_or(-1) + 1);
        *stream_it = separator;
    }
    else
    {
        stream_it.write_number(get<1>(mate));
        *stream_it = separator;
    }

    stream_it.write_number(get<2>(mate));
    *stream_it = separator;

    write_range_or_asterisk(stream_it, seq);
    *stream_it = separator;

    write_range_or_asterisk(stream_it, qual);

    write_tag_fields(stream_it, tag_dict, separator);

    stream_it.write_end_of_line(options.add_carriage_return);
}

/*!\brief Reads a list of values separated by comma as it is the case for SAM tag arrays.
 * \tparam value_type       The type of values to be stored in the tag array.
 *
 * \param[in, out] variant      A std::variant object to store the tag arrays.
 * \param[in, out] str          The string_view to parse.
 * \param[in]      value        A temporary value that determines the underlying type of the tag array.
 *
 * \details
 *
 * Reading the tags is done according to the official
 * [SAM format specifications](https://samtools.github.io/hts-specs/SAMv1.pdf).
 *
 * The function throws a seqan3::format_error if any unknown tag type was encountered. It will also fail if the
 * format is not in a correct state (e.g. required fields are not given), but throwing might occur downstream of
 * the actual error.
 */
template <arithmetic value_type>
inline void format_sam::read_sam_dict_vector(seqan3::detail::sam_tag_variant & variant,
                                             std::string_view const str,
                                             value_type value)
{
    std::vector<value_type> tmp_vector{};
    size_t start_pos{0};
    size_t end_pos{0};

    while (start_pos != std::string_view::npos)
    {
        end_pos = str.find(',', start_pos);
        auto end = (end_pos == std::string_view::npos) ? str.end() : str.begin() + end_pos;
        read_arithmetic_field(std::string_view{str.begin() + start_pos, end}, value);
        tmp_vector.push_back(value);

        start_pos = (end_pos == std::string_view::npos) ? end_pos : end_pos + 1;
    }
    variant = std::move(tmp_vector);
}

/*!\brief Reads a list of byte pairs as it is the case for SAM tag byte arrays.
 * \param[in, out] variant      A std::variant object to store the tag arrays.
 * \param[in, out] str          The string_view to parse.
 *
 * \details
 *
 * Reading the byte tags is done according to the official
 * [SAM format specifications](https://samtools.github.io/hts-specs/SAMv1.pdf).
 *
 * The function throws a seqan3::format_error if there was an uneven number of bytes.
 */
inline void format_sam::read_sam_byte_vector(seqan3::detail::sam_tag_variant & variant, std::string_view const str)
{
    std::vector<std::byte> tmp_vector{};
    // std::from_chars cannot directly parse into a std::byte
    uint8_t dummy_byte{};

    if (str.size() % 2 != 0)
        throw format_error{"[CORRUPTED SAM FILE]  Hexadecimal tag must have even number of digits."};

    // H encodes bytes in a hexadecimal format. Two hex values are stored for each byte as characters.
    // E.g., '1' and 'A' need one byte each and are read as `\x1A`, which is 27 in decimal.
    for (auto hex_begin = str.begin(), hex_end = str.begin() + 2; hex_begin != str.end(); hex_begin += 2, hex_end += 2)
    {
        auto res = std::from_chars(hex_begin, hex_end, dummy_byte, 16);

        if (res.ec == std::errc::invalid_argument)
            throw format_error{std::string("[CORRUPTED SAM FILE] The string '") + std::string(hex_begin, hex_end)
                               + "' could not be cast into type uint8_t."};

        if (res.ec == std::errc::result_out_of_range)
            throw format_error{std::string("[CORRUPTED SAM FILE] Casting '") + std::string(str)
                               + "' into type uint8_t would cause an overflow."};

        tmp_vector.push_back(std::byte{dummy_byte});
    }

    variant = std::move(tmp_vector);
}

/*!\brief Reads the optional tag fields into the seqan3::sam_tag_dictionary.
 * \param[in, out] tag_str      The string_view to parse for the sam_tag_dictionary entries.
 * \param[in, out] target       The seqan3::sam_tag_dictionary to store the tag information.
 *
 * \throws seqan3::format_error if any unexpected character or format is encountered.
 *
 * \details
 *
 * Reading the tags is done according to the official
 * [SAM format specifications](https://samtools.github.io/hts-specs/SAMv1.pdf).
 *
 * The function throws a seqan3::format_error if any unknown tag type was encountered. It will also fail if the
 * format is not in a correct state (e.g. required fields are not given), but throwing might occur downstream of
 * the actual error.
 */
inline void format_sam::read_sam_dict(std::string_view const tag_str, sam_tag_dictionary & target)
{
    /* Every SAM tag has the format "[TAG]:[TYPE_ID]:[VALUE]", where TAG is a two letter
       name tag which is converted to a unique integer identifier and TYPE_ID is one character in [A,i,Z,H,B,f]
       describing the type for the upcoming VALUES. If TYPE_ID=='B' it signals an array of comma separated
       VALUE's and the inner value type is identified by the character following ':', one of [cCsSiIf].
    */
    assert(tag_str.size() > 5);

    uint16_t tag = static_cast<uint16_t>(tag_str[0]) << 8;
    tag += static_cast<uint16_t>(tag_str[1]);

    char type_id = tag_str[3];

    switch (type_id)
    {
    case 'A': // char
    {
        assert(tag_str.size() == 6);
        target[tag] = tag_str[5];
        break;
    }
    case 'i': // int32_t
    {
        int32_t tmp;
        read_arithmetic_field(tag_str.substr(5), tmp);
        target[tag] = tmp;
        break;
    }
    case 'f': // float
    {
        float tmp;
        read_arithmetic_field(tag_str.substr(5), tmp);
        target[tag] = tmp;
        break;
    }
    case 'Z': // string
    {
        target[tag] = std::string{tag_str.substr(5)};
        break;
    }
    case 'H':
    {
        read_sam_byte_vector(target[tag], tag_str.substr(5));
        break;
    }
    case 'B': // Array. Value type depends on second char [cCsSiIf]
    {
        assert(tag_str.size() > 6);
        char array_value_type_id = tag_str[5];

        switch (array_value_type_id)
        {
        case 'c': // int8_t
            read_sam_dict_vector(target[tag], tag_str.substr(7), int8_t{});
            break;
        case 'C': // uint8_t
            read_sam_dict_vector(target[tag], tag_str.substr(7), uint8_t{});
            break;
        case 's': // int16_t
            read_sam_dict_vector(target[tag], tag_str.substr(7), int16_t{});
            break;
        case 'S': // uint16_t
            read_sam_dict_vector(target[tag], tag_str.substr(7), uint16_t{});
            break;
        case 'i': // int32_t
            read_sam_dict_vector(target[tag], tag_str.substr(7), int32_t{});
            break;
        case 'I': // uint32_t
            read_sam_dict_vector(target[tag], tag_str.substr(7), uint32_t{});
            break;
        case 'f': // float
            read_sam_dict_vector(target[tag], tag_str.substr(7), float{});
            break;
        default:
            throw format_error{std::string("The first character in the numerical ")
                               + "id of a SAM tag must be one of [cCsSiIf] but '" + array_value_type_id
                               + "' was given."};
        }
        break;
    }
    default:
        throw format_error{std::string("The second character in the numerical id of a "
                                       "SAM tag ([TAG]:[TYPE_ID]:[VALUE]) must be one of [A,i,Z,H,B,f] but '")
                           + type_id + "' was given."};
    }
}

/*!\brief Writes a field value to the stream.
 * \tparam stream_it_t The stream iterator type.
 * \tparam field_type  The type of the field value. Must model std::ranges::forward_range.
 *
 * \param[in,out] stream_it   The stream iterator to print to.
 * \param[in]     field_value The value to print.
 */
template <typename stream_it_t, std::ranges::forward_range field_type>
inline void format_sam::write_range_or_asterisk(stream_it_t & stream_it, field_type && field_value)
{
    if (std::ranges::empty(field_value))
    {
        *stream_it = '*';
    }
    else
    {
        if constexpr (std::same_as<std::remove_cvref_t<std::ranges::range_reference_t<field_type>>, char>)
            stream_it.write_range(field_value);
        else // convert from alphabets to their character representation
            stream_it.write_range(field_value | views::to_char);
    }
}

/*!\brief Writes a field value to the stream.
 * \tparam stream_it_t The stream iterator type.
 *
 * \param[in,out] stream_it   The stream iterator to print to.
 * \param[in]     field_value The value to print; a null-terminated CString.
 */
template <typename stream_it_t>
inline void format_sam::write_range_or_asterisk(stream_it_t & stream_it, char const * const field_value)
{
    write_range_or_asterisk(stream_it, std::string_view{field_value});
}

/*!\brief Writes the optional fields of the seqan3::sam_tag_dictionary.
 * \tparam stream_it_t      The stream iterator's type.
 *
 * \param[in,out] stream_it The stream iterator to print to.
 * \param[in]     tag_dict  The tag dictionary to print.
 * \param[in]     separator The field separator to append.
 */
template <typename stream_it_t>
inline void
format_sam::write_tag_fields(stream_it_t & stream_it, sam_tag_dictionary const & tag_dict, char const separator)
{
    auto const stream_variant_fn = [&stream_it](auto && arg) // helper to print a std::variant
    {
        using T = std::remove_cvref_t<decltype(arg)>;

        if constexpr (std::ranges::input_range<T>)
        {
            if constexpr (std::same_as<std::remove_cvref_t<std::ranges::range_reference_t<T>>, char>)
            {
                stream_it.write_range(arg);
            }
            else if constexpr (std::same_as<std::remove_cvref_t<std::ranges::range_reference_t<T>>, std::byte>)
            {
                if (!std::ranges::empty(arg))
                {
                    stream_it.write_number(std::to_integer<uint8_t>(*std::ranges::begin(arg)));

                    for (auto && elem : arg | std::views::drop(1))
                    {
                        *stream_it = ',';
                        stream_it.write_number(std::to_integer<uint8_t>(elem));
                    }
                }
            }
            else
            {
                if (!std::ranges::empty(arg))
                {
                    stream_it.write_number(*std::ranges::begin(arg));

                    for (auto && elem : arg | std::views::drop(1))
                    {
                        *stream_it = ',';
                        stream_it.write_number(elem);
                    }
                }
            }
        }
        else if constexpr (std::same_as<std::remove_cvref_t<T>, char>)
        {
            *stream_it = arg;
        }
        else // number
        {
            stream_it.write_number(arg);
        }
    };

    for (auto & [tag, variant] : tag_dict)
    {
        *stream_it = separator;

        char const char0 = tag / 256;
        char const char1 = tag % 256;

        *stream_it = char0;
        *stream_it = char1;
        *stream_it = ':';
        *stream_it = detail::sam_tag_type_char[variant.index()];
        *stream_it = ':';

        if (detail::sam_tag_type_char_extra[variant.index()] != '\0')
        {
            *stream_it = detail::sam_tag_type_char_extra[variant.index()];
            *stream_it = ',';
        }

        std::visit(stream_variant_fn, variant);
    }
}

} // namespace seqan3
