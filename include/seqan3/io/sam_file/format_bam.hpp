// -----------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2021, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2021, Knut Reinert & MPI für molekulare Genetik
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at: https://github.com/seqan/seqan3/blob/master/LICENSE.md
// -----------------------------------------------------------------------------------------------------

/*!\file
 * \brief Provides the seqan3::format_bam.
 * \author Svenja Mehringer <svenja.mehringer AT fu-berlin.de>
 */

#pragma once

#include <seqan3/std/bit>
#include <iterator>
#include <seqan3/std/ranges>
#include <string>
#include <vector>

#include <seqan3/alphabet/nucleotide/dna16sam.hpp>
#include <seqan3/core/debug_stream/optional.hpp>
#include <seqan3/io/sam_file/detail/cigar.hpp>
#include <seqan3/io/sam_file/detail/format_sam_base.hpp>
#include <seqan3/io/sam_file/header.hpp>
#include <seqan3/io/sam_file/input_options.hpp>
#include <seqan3/io/sam_file/sam_flag.hpp>
#include <seqan3/io/sam_file/sam_tag_dictionary.hpp>
#include <seqan3/io/stream/detail/fast_ostreambuf_iterator.hpp>
#include <seqan3/io/views/detail/istreambuf_view.hpp>
#include <seqan3/io/views/detail/take_exactly_view.hpp>
#include <seqan3/utility/views/slice.hpp>

namespace seqan3
{

/*!\brief       The BAM format.
 * \implements  AlignmentFileFormat
 * \ingroup io_sam_file
 *
 * \details
 *
 * The BAM format is the binary version of the SAM format:
 *
 * \copydetails seqan3::format_sam
 */
class format_bam : private detail::format_sam_base
{
public:
    /*!\name Constructors, destructor and assignment
     * \{
     */
    // string_buffer is of type std::string and has some problems with pre-C++11 ABI
    format_bam() = default; //!< Defaulted.
    format_bam(format_bam const &) = default; //!< Defaulted.
    format_bam & operator=(format_bam const &) = default; //!< Defaulted.
    format_bam(format_bam &&) = default; //!< Defaulted.
    format_bam & operator=(format_bam &&) = default; //!< Defaulted.
    ~format_bam() = default; //!< Defaulted.

    //!\}

    //!\brief The valid file extensions for this format; note that you can modify this value.
    static inline std::vector<std::string> file_extensions
    {
        { "bam" }
    };

protected:
    template <typename stream_type,     // constraints checked by file
              typename seq_legal_alph_type,
              typename ref_seqs_type,
              typename ref_ids_type,
              typename stream_pos_type,
              typename seq_type,
              typename id_type,
              typename offset_type,
              typename ref_seq_type,
              typename ref_id_type,
              typename ref_offset_type,
              typename align_type,
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
                               offset_type & offset,
                               ref_seq_type & SEQAN3_DOXYGEN_ONLY(ref_seq),
                               ref_id_type & ref_id,
                               ref_offset_type & ref_offset,
                               align_type & align,
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
              typename align_type,
              typename cigar_type,
              typename qual_type,
              typename mate_type,
              typename tag_dict_type>
    void write_alignment_record([[maybe_unused]] stream_type &  stream,
                                [[maybe_unused]] sam_file_output_options const & options,
                                [[maybe_unused]] header_type && header,
                                [[maybe_unused]] seq_type && seq,
                                [[maybe_unused]] qual_type && qual,
                                [[maybe_unused]] id_type && id,
                                [[maybe_unused]] int32_t const offset,
                                [[maybe_unused]] ref_seq_type && SEQAN3_DOXYGEN_ONLY(ref_seq),
                                [[maybe_unused]] ref_id_type && ref_id,
                                [[maybe_unused]] std::optional<int32_t> ref_offset,
                                [[maybe_unused]] align_type && align,
                                [[maybe_unused]] cigar_type && cigar_vector,
                                [[maybe_unused]] sam_flag const flag,
                                [[maybe_unused]] uint8_t const mapq,
                                [[maybe_unused]] mate_type && mate,
                                [[maybe_unused]] tag_dict_type && tag_dict,
                                [[maybe_unused]] double SEQAN3_DOXYGEN_ONLY(e_value),
                                [[maybe_unused]] double SEQAN3_DOXYGEN_ONLY(bit_score));

private:
    //!\brief A variable that tracks whether the content of header has been read or not.
    bool header_was_read{false};

    //!\brief Local buffer to read into while avoiding reallocation.
    std::string string_buffer{};

    //!\brief Stores all fixed length variables which can be read/written directly by reinterpreting the binary stream.
    struct alignment_record_core
    {   // naming corresponds to official SAM/BAM specifications
        int32_t block_size;     //!< The size in bytes of the whole BAM record.
        int32_t refID;          //!< The reference id the read was mapped to.
        int32_t pos;            //!< The begin position of the alignment.
        uint32_t l_read_name:8; //!< The length of the read name including the \0 character.
        uint32_t mapq:8;        //!< The mapping quality.
        uint32_t bin:16;        //!< The bin number.
        uint32_t n_cigar_op:16; //!< The number of cigar operations of the alignment.
        sam_flag flag;          //!< The flag value (uint16_t enum).
        int32_t l_seq;          //!< The number of bases of the read sequence.
        int32_t next_refID;     //!< The reference id of the mate.
        int32_t next_pos;       //!< The begin position of the mate alignment.
        int32_t tlen;           //!< The template length of the read and its mate.
    };

    //!\brief Converts a cigar op character to the rank according to the official BAM specifications.
    static constexpr std::array<uint8_t, 256> char_to_sam_rank
    {
        [] () constexpr
        {
            std::array<uint8_t, 256> ret{};

            using index_t = std::make_unsigned_t<char>;

            // ret['M'] = 0; set anyway by initialization
            ret[static_cast<index_t>('I')] = 1;
            ret[static_cast<index_t>('D')] = 2;
            ret[static_cast<index_t>('N')] = 3;
            ret[static_cast<index_t>('S')] = 4;
            ret[static_cast<index_t>('H')] = 5;
            ret[static_cast<index_t>('P')] = 6;
            ret[static_cast<index_t>('=')] = 7;
            ret[static_cast<index_t>('X')] = 8;

            return ret;
        }()
    };

    //!\brief Computes the bin number for a given region [beg, end), copied from the official SAM specifications.
    static uint16_t reg2bin(int32_t beg, int32_t end) noexcept
    {
        --end;
        if (beg >> 14 == end >> 14) return ((1 << 15) - 1) / 7 + (beg >> 14);
        if (beg >> 17 == end >> 17) return ((1 << 12) - 1) / 7 + (beg >> 17);
        if (beg >> 20 == end >> 20) return ((1 << 9)  - 1) / 7 + (beg >> 20);
        if (beg >> 23 == end >> 23) return ((1 << 6)  - 1) / 7 + (beg >> 23);
        if (beg >> 26 == end >> 26) return ((1 << 3)  - 1) / 7 + (beg >> 26);
        return 0;
    }

    using format_sam_base::read_field; // inherit read_field functions from format_base explicitly

    /*!\brief Reads a arithmetic field from binary stream by directly reinterpreting the bits.
     * \tparam stream_view_type  The type of the stream as a view.
     * \tparam number_type       The type of number to parse; must model std::integral.
     * \param[in, out] stream_view  The stream view to read from.
     * \param[out]     target       An integral value to store the parsed value in.
     */
    template <typename stream_view_type, std::integral number_type>
    void read_field(stream_view_type && stream_view, number_type & target)
    {
        std::ranges::copy_n(std::ranges::begin(stream_view), sizeof(target), reinterpret_cast<char *>(&target));
    }

    /*!\brief Reads a float field from binary stream by directly reinterpreting the bits.
     * \tparam stream_view_type  The type of the stream as a view.
     * \param[in, out] stream_view  The stream view to read from.
     * \param[out]     target       An float value to store the parsed value in.
     */
    template <typename stream_view_type>
    void read_field(stream_view_type && stream_view, float & target)
    {
        std::ranges::copy_n(std::ranges::begin(stream_view), sizeof(int32_t), reinterpret_cast<char *>(&target));
    }

    /*!\brief Delegate parsing of std::optional types to parsing of the inner value type.
     * \tparam stream_view_type     The type of the stream as a view.
     * \tparam optional_value_type  The inner type of a the std::optional type of \p target.
     *
     * \param[in, out] stream_view  The stream view to iterate over.
     * \param[out]     target       The std::optional object to store the parsed value.
     *
     * \throws seqan3::format_error if the character sequence in stream_view cannot be successfully converted to a value
     *         of type target_type.
     */
    template <typename stream_view_type, typename optional_value_type>
    void read_field(stream_view_type && stream_view, std::optional<optional_value_type> & target)
    {
        optional_value_type tmp;
        read_field(std::forward<stream_view_type>(stream_view), tmp);
        target = tmp;
    }

    template <typename stream_view_type, typename value_type>
    void read_sam_dict_vector(seqan3::detail::sam_tag_variant & variant,
                              stream_view_type && stream_view,
                              value_type const & SEQAN3_DOXYGEN_ONLY(value));

    template <typename stream_view_type>
    void read_field(stream_view_type && stream_view, sam_tag_dictionary & target);

    template <typename cigar_input_type>
    auto parse_binary_cigar(cigar_input_type && cigar_input, uint16_t n_cigar_op) const;

    static std::string get_tag_dict_str(sam_tag_dictionary const & tag_dict);
};

//!\copydoc sam_file_input_format::read_alignment_record
template <typename stream_type,     // constraints checked by file
          typename seq_legal_alph_type,
          typename ref_seqs_type,
          typename ref_ids_type,
          typename stream_pos_type,
          typename seq_type,
          typename id_type,
          typename offset_type,
          typename ref_seq_type,
          typename ref_id_type,
          typename ref_offset_type,
          typename align_type,
          typename cigar_type,
          typename flag_type,
          typename mapq_type,
          typename qual_type,
          typename mate_type,
          typename tag_dict_type,
          typename e_value_type,
          typename bit_score_type>
inline void format_bam::read_alignment_record(stream_type & stream,
                                              sam_file_input_options<seq_legal_alph_type> const & SEQAN3_DOXYGEN_ONLY(options),
                                              ref_seqs_type & ref_seqs,
                                              sam_file_header<ref_ids_type> & header,
                                              stream_pos_type & position_buffer,
                                              seq_type & seq,
                                              qual_type & qual,
                                              id_type & id,
                                              offset_type & offset,
                                              ref_seq_type & SEQAN3_DOXYGEN_ONLY(ref_seq),
                                              ref_id_type & ref_id,
                                              ref_offset_type & ref_offset,
                                              align_type & align,
                                              cigar_type & cigar_vector,
                                              flag_type & flag,
                                              mapq_type & mapq,
                                              mate_type & mate,
                                              tag_dict_type & tag_dict,
                                              e_value_type & SEQAN3_DOXYGEN_ONLY(e_value),
                                              bit_score_type & SEQAN3_DOXYGEN_ONLY(bit_score))
{
    static_assert(detail::decays_to_ignore_v<ref_offset_type> ||
                  detail::is_type_specialisation_of_v<ref_offset_type, std::optional>,
                  "The ref_offset must be a specialisation of std::optional.");

    static_assert(detail::decays_to_ignore_v<mapq_type> || std::same_as<mapq_type, uint8_t>,
                  "The type of field::mapq must be uint8_t.");

    static_assert(detail::decays_to_ignore_v<flag_type> || std::same_as<flag_type, sam_flag>,
                  "The type of field::flag must be seqan3::sam_flag.");

    auto stream_view = seqan3::detail::istreambuf(stream);

    // these variables need to be stored to compute the ALIGNMENT
    [[maybe_unused]] int32_t offset_tmp{};
    [[maybe_unused]] int32_t soft_clipping_end{};
    [[maybe_unused]] std::vector<cigar> tmp_cigar_vector{};
    [[maybe_unused]] int32_t ref_length{0}, seq_length{0}; // length of aligned part for ref and query

    // Header
    // -------------------------------------------------------------------------------------------------------------
    if (!header_was_read)
    {
        // magic BAM string
        if (!std::ranges::equal(stream_view | detail::take_exactly_or_throw(4), std::string_view{"BAM\1"}))
            throw format_error{"File is not in BAM format."};

        int32_t l_text{}; // length of header text including \0 character
        int32_t n_ref{}; // number of reference sequences
        int32_t l_name{}; // 1 + length of reference name including \0 character
        int32_t l_ref{}; // length of reference sequence

        read_field(stream_view, l_text);

        if (l_text > 0) // header text is present
            read_header(stream_view | detail::take_exactly_or_throw(l_text), header, ref_seqs);

        read_field(stream_view, n_ref);

        for (int32_t ref_idx = 0; ref_idx < n_ref; ++ref_idx)
        {
            read_field(stream_view, l_name);

            string_buffer.resize(l_name - 1);
            std::ranges::copy_n(std::ranges::begin(stream_view), l_name - 1, string_buffer.data()); // copy without \0 character
            std::ranges::next(std::ranges::begin(stream_view)); // skip \0 character

            read_field(stream_view, l_ref);

            if constexpr (detail::decays_to_ignore_v<ref_seqs_type>) // no reference information given
            {
                // If there was no header text, we parse reference sequences block as header information
                if (l_text == 0)
                {
                    auto & reference_ids = header.ref_ids();
                    // put the length of the reference sequence into ref_id_info
                    header.ref_id_info.emplace_back(l_ref, "");
                    // put the reference name into reference_ids
                    reference_ids.push_back(string_buffer);
                    // assign the reference name an ascending reference id (starts at index 0).
                    header.ref_dict.emplace(reference_ids.back(), reference_ids.size() - 1);
                    continue;
                }
            }

            auto id_it = header.ref_dict.find(string_buffer);

            // sanity checks of reference information to existing header object:
            if (id_it == header.ref_dict.end()) // [unlikely]
            {
                throw format_error{detail::to_string("Unknown reference name '" + string_buffer +
                                                     "' found in BAM file header (header.ref_ids():",
                                                     header.ref_ids(), ").")};
            }
            else if (id_it->second != ref_idx) // [unlikely]
            {
                throw format_error{detail::to_string("Reference id '", string_buffer, "' at position ", ref_idx,
                                                     " does not correspond to the position ", id_it->second,
                                                     " in the header (header.ref_ids():", header.ref_ids(), ").")};
            }
            else if (std::get<0>(header.ref_id_info[id_it->second]) != l_ref) // [unlikely]
            {
                throw format_error{"Provided reference has unequal length as specified in the header."};
            }
        }

        header_was_read = true;

        if (std::ranges::begin(stream_view) == std::ranges::end(stream_view)) // no records follow
            return;
    }

    // read alignment record into buffer
    // -------------------------------------------------------------------------------------------------------------
    position_buffer = stream.tellg();
    alignment_record_core core;
    std::ranges::copy(stream_view | detail::take_exactly_or_throw(sizeof(core)), reinterpret_cast<char *>(&core));

    if (core.refID >= static_cast<int32_t>(header.ref_ids().size()) || core.refID < -1) // [[unlikely]]
    {
        throw format_error{detail::to_string("Reference id index '", core.refID, "' is not in range of ",
                                             "header.ref_ids(), which has size ", header.ref_ids().size(), ".")};
    }
    else if (core.refID > -1) // not unmapped
    {
        ref_id = core.refID;                                                   // field::ref_id
    }

    flag = core.flag;                                                          // field::flag
    mapq = core.mapq;                                                          // field::mapq

    if (core.pos > -1) // [[likely]]
        ref_offset = core.pos;                                                 // field::ref_offset

    if constexpr (!detail::decays_to_ignore_v<mate_type>)                      // field::mate
    {
        if (core.next_refID > -1)
            get<0>(mate) = core.next_refID;

        if (core.next_pos > -1) // [[likely]]
            get<1>(mate) = core.next_pos;

        get<2>(mate) = core.tlen;
    }

    // read id
    // -------------------------------------------------------------------------------------------------------------
    read_field(stream_view | detail::take_exactly_or_throw(core.l_read_name - 1), id); // field::id
    std::ranges::next(std::ranges::begin(stream_view)); // skip '\0'

    // read cigar string
    // -------------------------------------------------------------------------------------------------------------
    if constexpr (!detail::decays_to_ignore_v<align_type> || !detail::decays_to_ignore_v<cigar_type>)
    {
        std::tie(tmp_cigar_vector, ref_length, seq_length) = parse_binary_cigar(stream_view, core.n_cigar_op);
        transfer_soft_clipping_to(tmp_cigar_vector, offset_tmp, soft_clipping_end);
        // the actual cigar_vector is swapped with tmp_cigar_vector at the end to avoid copying
    }
    else
    {
        detail::consume(stream_view | detail::take_exactly_or_throw(core.n_cigar_op * 4));
    }

    offset = offset_tmp;

    // read sequence
    // -------------------------------------------------------------------------------------------------------------
    if (core.l_seq > 0) // sequence information is given
    {
        auto seq_stream = stream_view
                        | detail::take_exactly_or_throw(core.l_seq / 2) // one too short if uneven
                        | std::views::transform([] (char c) -> std::pair<dna16sam, dna16sam>
                          {
                              return {dna16sam{}.assign_rank(std::min(15, static_cast<uint8_t>(c) >> 4)),
                                      dna16sam{}.assign_rank(std::min(15, static_cast<uint8_t>(c) & 0x0f))};
                          });

        if constexpr (detail::decays_to_ignore_v<seq_type>)
        {
            auto skip_sequence_bytes = [&] ()
            {
                detail::consume(seq_stream);
                if (core.l_seq & 1)
                    std::ranges::next(std::ranges::begin(stream_view));
            };

            if constexpr (!detail::decays_to_ignore_v<align_type>)
            {
                static_assert(sequence_container<std::remove_reference_t<decltype(get<1>(align))>>,
                              "If you want to read ALIGNMENT but not SEQ, the alignment"
                              " object must store a sequence container at the second (query) position.");

                if (!tmp_cigar_vector.empty()) // only parse alignment if cigar information was given
                {
                    assert(core.l_seq == (seq_length + offset_tmp + soft_clipping_end)); // sanity check
                    using alph_t = std::ranges::range_value_t<decltype(get<1>(align))>;
                    constexpr auto from_dna16 = detail::convert_through_char_representation<alph_t, dna16sam>;

                    get<1>(align).reserve(seq_length);

                    auto tmp_iter = std::ranges::begin(seq_stream);
                    std::ranges::advance(tmp_iter, offset_tmp / 2); // skip soft clipped bases at the beginning

                    if (offset_tmp & 1)
                    {
                        get<1>(align).push_back(from_dna16[to_rank(get<1>(*tmp_iter))]);
                        ++tmp_iter;
                        --seq_length;
                    }

                    for (size_t i = (seq_length / 2); i > 0; --i)
                    {
                        get<1>(align).push_back(from_dna16[to_rank(get<0>(*tmp_iter))]);
                        get<1>(align).push_back(from_dna16[to_rank(get<1>(*tmp_iter))]);
                        ++tmp_iter;
                    }

                    if (seq_length & 1)
                    {
                        get<1>(align).push_back(from_dna16[to_rank(get<0>(*tmp_iter))]);
                        ++tmp_iter;
                        --soft_clipping_end;
                    }

                    std::ranges::advance(tmp_iter, (soft_clipping_end + !(seq_length & 1)) / 2);
                }
                else
                {
                    skip_sequence_bytes();
                    get<1>(align) = std::remove_reference_t<decltype(get<1>(align))>{}; // assign empty container
                }
            }
            else
            {
                skip_sequence_bytes();
            }
        }
        else
        {
            using alph_t = std::ranges::range_value_t<decltype(seq)>;
            constexpr auto from_dna16 = detail::convert_through_char_representation<alph_t, dna16sam>;

            for (auto [d1, d2] : seq_stream)
            {
                seq.push_back(from_dna16[to_rank(d1)]);
                seq.push_back(from_dna16[to_rank(d2)]);
            }

            if (core.l_seq & 1)
            {
                dna16sam d = dna16sam{}.assign_rank(std::min(15, static_cast<uint8_t>(*std::ranges::begin(stream_view)) >> 4));
                seq.push_back(from_dna16[to_rank(d)]);
                std::ranges::next(std::ranges::begin(stream_view));
            }

            if constexpr (!detail::decays_to_ignore_v<align_type>)
            {
                assign_unaligned(get<1>(align),
                                 seq | views::slice(static_cast<std::ranges::range_difference_t<seq_type>>(offset_tmp),
                                                    std::ranges::distance(seq) - soft_clipping_end));
            }
        }
    }

    // read qual string
    // -------------------------------------------------------------------------------------------------------------
    read_field(stream_view | detail::take_exactly_or_throw(core.l_seq)
                           | std::views::transform([] (char chr) { return static_cast<char>(chr + 33); }), qual);

    // All remaining optional fields if any: SAM tags dictionary
    // -------------------------------------------------------------------------------------------------------------
    int32_t remaining_bytes = core.block_size - (sizeof(alignment_record_core) - 4/*block_size excluded*/) -
                              core.l_read_name - core.n_cigar_op * 4 - (core.l_seq + 1) / 2 - core.l_seq;
    assert(remaining_bytes >= 0);
    auto tags_view = stream_view | detail::take_exactly_or_throw(remaining_bytes);

    while (tags_view.size() > 0)
        read_field(tags_view, tag_dict);

    // DONE READING - wrap up
    // -------------------------------------------------------------------------------------------------------------
    if constexpr (!detail::decays_to_ignore_v<align_type> || !detail::decays_to_ignore_v<cigar_type>)
    {
        // Check cigar, if it matches ‘kSmN’, where ‘k’ equals lseq, ‘m’ is the reference sequence length in the
        // alignment, and ‘S’ and ‘N’ are the soft-clipping and reference-clip, then the cigar string was larger
        // than 65535 operations and is stored in the sam_tag_dictionary (tag GC).
        if (core.l_seq != 0 && offset_tmp == core.l_seq)
        {
            if constexpr (detail::decays_to_ignore_v<tag_dict_type> | detail::decays_to_ignore_v<seq_type>)
            { // maybe only throw in debug mode and otherwise return an empty alignment?
                throw format_error{detail::to_string("The cigar string '", offset_tmp, "S", ref_length,
                                   "N' suggests that the cigar string exceeded 65535 elements and was therefore ",
                                   "stored in the optional field CG. You need to read in the field::tags and "
                                   "field::seq in order to access this information.")};
            }
            else
            {
                auto it = tag_dict.find("CG"_tag);

                if (it == tag_dict.end())
                    throw format_error{detail::to_string("The cigar string '", offset_tmp, "S", ref_length,
                                   "N' suggests that the cigar string exceeded 65535 elements and was therefore ",
                                   "stored in the optional field CG but this tag is not present in the given ",
                                   "record.")};

                auto cigar_view = std::views::all(std::get<std::string>(it->second));
                std::tie(tmp_cigar_vector, ref_length, seq_length) = detail::parse_cigar(cigar_view);
                offset_tmp = soft_clipping_end = 0;
                transfer_soft_clipping_to(tmp_cigar_vector, offset_tmp, soft_clipping_end);
                tag_dict.erase(it); // remove redundant information

                if constexpr (!detail::decays_to_ignore_v<align_type>)
                {
                    assign_unaligned(get<1>(align),
                                     seq | views::slice(static_cast<std::ranges::range_difference_t<seq_type>>(offset_tmp),
                                                        std::ranges::distance(seq) - soft_clipping_end));
                }
            }
        }
    }

    // Alignment object construction
    if constexpr (!detail::decays_to_ignore_v<align_type>)
        construct_alignment(align, tmp_cigar_vector, core.refID, ref_seqs, core.pos, ref_length); // inherited from SAM

    if constexpr (!detail::decays_to_ignore_v<cigar_type>)
        std::swap(cigar_vector, tmp_cigar_vector);
}

//!\copydoc sam_file_output_format::write_alignment_record
template <typename stream_type,
          typename header_type,
          typename seq_type,
          typename id_type,
          typename ref_seq_type,
          typename ref_id_type,
          typename align_type,
          typename cigar_type,
          typename qual_type,
          typename mate_type,
          typename tag_dict_type>
inline void format_bam::write_alignment_record([[maybe_unused]] stream_type &  stream,
                                               [[maybe_unused]] sam_file_output_options const & options,
                                               [[maybe_unused]] header_type && header,
                                               [[maybe_unused]] seq_type && seq,
                                               [[maybe_unused]] qual_type && qual,
                                               [[maybe_unused]] id_type && id,
                                               [[maybe_unused]] int32_t const offset,
                                               [[maybe_unused]] ref_seq_type && SEQAN3_DOXYGEN_ONLY(ref_seq),
                                               [[maybe_unused]] ref_id_type && ref_id,
                                               [[maybe_unused]] std::optional<int32_t> ref_offset,
                                               [[maybe_unused]] align_type && align,
                                               [[maybe_unused]] cigar_type && cigar_vector,
                                               [[maybe_unused]] sam_flag const flag,
                                               [[maybe_unused]] uint8_t const mapq,
                                               [[maybe_unused]] mate_type && mate,
                                               [[maybe_unused]] tag_dict_type && tag_dict,
                                               [[maybe_unused]] double SEQAN3_DOXYGEN_ONLY(e_value),
                                               [[maybe_unused]] double SEQAN3_DOXYGEN_ONLY(bit_score))
{
    // ---------------------------------------------------------------------
    // Type Requirements (as static asserts for user friendliness)
    // ---------------------------------------------------------------------
    static_assert((std::ranges::forward_range<seq_type>        &&
                  alphabet<std::ranges::range_reference_t<seq_type>>),
                  "The seq object must be a std::ranges::forward_range over "
                  "letters that model seqan3::alphabet.");

    static_assert((std::ranges::forward_range<id_type>         &&
                  alphabet<std::ranges::range_reference_t<id_type>>),
                  "The id object must be a std::ranges::forward_range over "
                  "letters that model seqan3::alphabet.");

    static_assert((std::ranges::forward_range<ref_seq_type>    &&
                  alphabet<std::ranges::range_reference_t<ref_seq_type>>),
                  "The ref_seq object must be a std::ranges::forward_range "
                  "over letters that model seqan3::alphabet.");

    if constexpr (!detail::decays_to_ignore_v<ref_id_type>)
    {
        static_assert((std::ranges::forward_range<ref_id_type> ||
                       std::integral<std::remove_reference_t<ref_id_type>> ||
                       detail::is_type_specialisation_of_v<std::remove_cvref_t<ref_id_type>, std::optional>),
                      "The ref_id object must be a std::ranges::forward_range "
                      "over letters that model seqan3::alphabet or an integral or a std::optional<integral>.");
    }

    static_assert(tuple_like<std::remove_cvref_t<align_type>>,
                  "The align object must be a std::pair of two ranges whose "
                  "value_type is comparable to seqan3::gap");

    static_assert((std::tuple_size_v<std::remove_cvref_t<align_type>> == 2 &&
                   std::equality_comparable_with<gap, std::ranges::range_reference_t<decltype(std::get<0>(align))>> &&
                   std::equality_comparable_with<gap, std::ranges::range_reference_t<decltype(std::get<1>(align))>>),
                  "The align object must be a std::pair of two ranges whose "
                  "value_type is comparable to seqan3::gap");

    static_assert((std::ranges::forward_range<qual_type>       &&
                   alphabet<std::ranges::range_reference_t<qual_type>>),
                  "The qual object must be a std::ranges::forward_range "
                  "over letters that model seqan3::alphabet.");

    static_assert(tuple_like<std::remove_cvref_t<mate_type>>,
                  "The mate object must be a std::tuple of size 3 with "
                  "1) a std::ranges::forward_range with a value_type modelling seqan3::alphabet, "
                  "2) a std::integral or std::optional<std::integral>, and "
                  "3) a std::integral.");

    static_assert(((std::ranges::forward_range<decltype(std::get<0>(mate))>     ||
                    std::integral<std::remove_cvref_t<decltype(std::get<0>(mate))>> ||
                    detail::is_type_specialisation_of_v<std::remove_cvref_t<decltype(std::get<0>(mate))>, std::optional>) &&
                  (std::integral<std::remove_cvref_t<decltype(std::get<1>(mate))>> ||
                   detail::is_type_specialisation_of_v<std::remove_cvref_t<decltype(std::get<1>(mate))>, std::optional>) &&
                  std::integral<std::remove_cvref_t<decltype(std::get<2>(mate))>>),
                  "The mate object must be a std::tuple of size 3 with "
                  "1) a std::ranges::forward_range with a value_type modelling seqan3::alphabet, "
                  "2) a std::integral or std::optional<std::integral>, and "
                  "3) a std::integral.");

    static_assert(std::same_as<std::remove_cvref_t<tag_dict_type>, sam_tag_dictionary>,
                  "The tag_dict object must be of type seqan3::sam_tag_dictionary.");

    if constexpr (detail::decays_to_ignore_v<header_type>)
    {
        throw format_error{"BAM can only be written with a header but you did not provide enough information! "
                           "You can either construct the output file with ref_ids and ref_seqs information and "
                           "the header will be created for you, or you can access the `header` member directly."};
    }
    else
    {
        // ---------------------------------------------------------------------
        // logical Requirements
        // ---------------------------------------------------------------------

        if (ref_offset.has_value() && (ref_offset.value() + 1) < 0)
            throw format_error{detail::to_string("The ref_offset object must be >= -1 but is: ", ref_offset)};

        detail::fast_ostreambuf_iterator stream_it{*stream.rdbuf()};

        // ---------------------------------------------------------------------
        // Writing the BAM Header on first call
        // ---------------------------------------------------------------------
        if (!header_was_written)
        {
            stream << "BAM\1";
            std::ostringstream os;
            write_header(os, options, header); // write SAM header to temporary stream to query the size.
            int32_t l_text{static_cast<int32_t>(os.str().size())};
            std::ranges::copy_n(reinterpret_cast<char *>(&l_text), 4, stream_it); // write read id

            stream  << os.str();

            int32_t n_ref{static_cast<int32_t>(header.ref_ids().size())};
            std::ranges::copy_n(reinterpret_cast<char *>(&n_ref), 4, stream_it); // write read id

            for (int32_t ridx = 0; ridx < n_ref; ++ridx)
            {
                int32_t l_name{static_cast<int32_t>(header.ref_ids()[ridx].size()) + 1}; // plus null character
                std::ranges::copy_n(reinterpret_cast<char *>(&l_name), 4, stream_it);    // write l_name
                // write reference name:
                std::ranges::copy(header.ref_ids()[ridx].begin(), header.ref_ids()[ridx].end(), stream_it);
                stream_it = '\0';
                // write reference sequence length:
                std::ranges::copy_n(reinterpret_cast<char *>(&get<0>(header.ref_id_info[ridx])), 4, stream_it);
            }

            header_was_written = true;
        }

        // ---------------------------------------------------------------------
        // Writing the Record
        // ---------------------------------------------------------------------
        int32_t ref_length{};

        // if alignment is non-empty, replace cigar_vector.
        // else, compute the ref_length from given cigar_vector which is needed to fill field `bin`.
        if (!std::ranges::empty(cigar_vector))
        {
            int32_t dummy_seq_length{};
            for (auto & [count, operation] : cigar_vector)
                detail::update_alignment_lengths(ref_length, dummy_seq_length, operation.to_char(), count);
        }
        else if (!std::ranges::empty(get<0>(align)) && !std::ranges::empty(get<1>(align)))
        {
            ref_length = std::ranges::distance(get<1>(align));

            // compute possible distance from alignment end to sequence end
            // which indicates soft clipping at the end.
            // This should be replaced by a free count_gaps function for
            // aligned sequences which is more efficient if possible.
            int32_t off_end{static_cast<int32_t>(std::ranges::distance(seq)) - offset};

            for (auto chr : get<1>(align))
                if (chr == gap{})
                    ++off_end;

            off_end -= ref_length;
            cigar_vector = detail::get_cigar_vector(align, offset, off_end);
        }

        if (cigar_vector.size() >= (1 << 16)) // must be written into the sam tag CG
        {
            tag_dict["CG"_tag] = detail::get_cigar_string(cigar_vector);
            cigar_vector.resize(2);
            cigar_vector[0] = cigar{static_cast<uint32_t>(std::ranges::distance(seq)), 'S'_cigar_operation};
            cigar_vector[1] = cigar{static_cast<uint32_t>(std::ranges::distance(get<1>(align))), 'N'_cigar_operation};
        }

        std::string tag_dict_binary_str = get_tag_dict_str(tag_dict);

        // Compute the value for the l_read_name field for the bam record.
        // This value is stored including a trailing `0`, so at most 254 characters of the id can be stored, since
        // the data type to store the value is uint8_t and 255 is the maximal size.
        // If the id is empty a '*' is written instead, i.e. the written id is never an empty string and stores at least
        // 2 bytes.
        uint8_t read_name_size = std::min<uint8_t>(std::ranges::distance(id), 254) + 1;
        read_name_size += static_cast<uint8_t>(read_name_size == 1); // need size two since empty id is stored as '*'.

        alignment_record_core core
        {
            /* block_size  */ 0,  // will be initialised right after
            /* refID       */ -1, // will be initialised right after
            /* pos         */ ref_offset.value_or(-1),
            /* l_read_name */ read_name_size,
            /* mapq        */ mapq,
            /* bin         */ reg2bin(ref_offset.value_or(-1), ref_length),
            /* n_cigar_op  */ static_cast<uint16_t>(cigar_vector.size()),
            /* flag        */ flag,
            /* l_seq       */ static_cast<int32_t>(std::ranges::distance(seq)),
            /* next_refId  */ -1, // will be initialised right after
            /* next_pos    */ get<1>(mate).value_or(-1),
            /* tlen        */ get<2>(mate)
        };

        auto check_and_assign_id_to = [&header] ([[maybe_unused]] auto & id_source,
                                                 [[maybe_unused]] auto & id_target)
        {
            using id_t = std::remove_reference_t<decltype(id_source)>;

            if constexpr (!detail::decays_to_ignore_v<id_t>)
            {
                if constexpr (std::integral<id_t>)
                {
                    id_target = id_source;
                }
                else if constexpr (detail::is_type_specialisation_of_v<id_t, std::optional>)
                {
                    id_target = id_source.value_or(-1);
                }
                else
                {
                    if (!std::ranges::empty(id_source)) // otherwise default will remain (-1)
                    {
                        auto id_it = header.ref_dict.end();

                        if constexpr (std::ranges::contiguous_range<decltype(id_source)> &&
                                      std::ranges::sized_range<decltype(id_source)> &&
                                      std::ranges::borrowed_range<decltype(id_source)>)
                        {
                            id_it = header.ref_dict.find(std::span{std::ranges::data(id_source),
                                                                   std::ranges::size(id_source)});
                        }
                        else
                        {
                            using header_ref_id_type = std::remove_reference_t<decltype(header.ref_ids()[0])>;

                            static_assert(implicitly_convertible_to<decltype(id_source), header_ref_id_type>,
                              "The ref_id type is not convertible to the reference id information stored in the "
                              "reference dictionary of the header object.");

                            id_it = header.ref_dict.find(id_source);
                        }

                        if (id_it == header.ref_dict.end())
                        {
                            throw format_error{detail::to_string("Unknown reference name '", id_source, "' could "
                                                                 "not be found in BAM header ref_dict: ",
                                                                 header.ref_dict, ".")};
                        }

                        id_target = id_it->second;
                    }
                }
            }
        };

        // initialise core.refID
        check_and_assign_id_to(ref_id, core.refID);

        // initialise core.next_refID
        check_and_assign_id_to(get<0>(mate), core.next_refID);

        // initialise core.block_size
        core.block_size = sizeof(core) - 4/*block_size excluded*/ +
                          core.l_read_name +
                          core.n_cigar_op * 4 +  // each int32_t has 4 bytes
                          (core.l_seq + 1) / 2 + // bitcompressed seq
                          core.l_seq +           // quality string
                          tag_dict_binary_str.size();

        std::ranges::copy_n(reinterpret_cast<char *>(&core), sizeof(core), stream_it);  // write core

        if (std::ranges::empty(id)) // empty id is represented as * for backward compatibility
            stream_it = '*';
        else
            std::ranges::copy_n(std::ranges::begin(id), core.l_read_name - 1, stream_it); // write read id
        stream_it = '\0';

        // write cigar
        for (auto [cigar_count, op] : cigar_vector)
        {
            cigar_count = cigar_count << 4;
            cigar_count |= static_cast<int32_t>(char_to_sam_rank[op.to_char()]);
            std::ranges::copy_n(reinterpret_cast<char *>(&cigar_count), 4, stream_it);
        }

        // write seq (bit-compressed: dna16sam characters go into one byte)
        using alph_t = std::ranges::range_value_t<seq_type>;
        constexpr auto to_dna16 = detail::convert_through_char_representation<dna16sam, alph_t>;

        auto sit = std::ranges::begin(seq);
        for (int32_t sidx = 0; sidx < ((core.l_seq & 1) ? core.l_seq - 1 : core.l_seq); ++sidx, ++sit)
        {
            uint8_t compressed_chr = to_rank(to_dna16[to_rank(*sit)]) << 4;
            ++sidx, ++sit;
            compressed_chr |= to_rank(to_dna16[to_rank(*sit)]);
            stream_it = static_cast<char>(compressed_chr);
        }

        if (core.l_seq & 1) // write one more
            stream_it = static_cast<char>(to_rank(to_dna16[to_rank(*sit)]) << 4);

        // write qual
        if (std::ranges::empty(qual))
        {
            auto v = views::repeat_n(static_cast<char>(255), core.l_seq);
            std::ranges::copy_n(v.begin(), core.l_seq, stream_it);
        }
        else
        {
            if (std::ranges::distance(qual) != core.l_seq)
                throw format_error{detail::to_string("Expected quality of same length as sequence with size ",
                                                     core.l_seq, ". Got quality with size ",
                                                     std::ranges::distance(qual), " instead.")};

            auto v = qual | std::views::transform([] (auto chr) { return static_cast<char>(to_rank(chr)); });
            std::ranges::copy_n(v.begin(), core.l_seq, stream_it);
        }

        // write optional fields
        stream << tag_dict_binary_str;
    } // if constexpr (!detail::decays_to_ignore_v<header_type>)
}

//!\copydoc seqan3::format_sam::read_sam_dict_vector
template <typename stream_view_type, typename value_type>
inline void format_bam::read_sam_dict_vector(seqan3::detail::sam_tag_variant & variant,
                                             stream_view_type && stream_view,
                                             value_type const & SEQAN3_DOXYGEN_ONLY(value))
{
    int32_t count;
    read_field(stream_view, count); // read length of vector
    std::vector<value_type> tmp_vector;
    tmp_vector.reserve(count);

    while (count > 0)
    {
        value_type tmp{};
        read_field(stream_view, tmp);
        tmp_vector.push_back(std::move(tmp));
        --count;
    }
    variant = std::move(tmp_vector);
}

/*!\brief Reads the optional tag fields into the seqan3::sam_tag_dictionary.
 * \tparam stream_view_type   The type of the stream as a view.
 *
 * \param[in, out] stream_view  The stream view to iterate over.
 * \param[out]     target       The seqan3::sam_tag_dictionary to store the tag information.
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
template <typename stream_view_type>
inline void format_bam::read_field(stream_view_type && stream_view, sam_tag_dictionary & target)
{
    /* Every BA< tag has the format "[TAG][TYPE_ID][VALUE]", where TAG is a two letter
       name tag which is converted to a unique integer identifier and TYPE_ID is one character in [A,i,Z,H,B,f]
       describing the type for the upcoming VALUES. If TYPE_ID=='B' it signals an array of
       VALUE's and the inner value type is identified by the next character, one of [cCsSiIf], followed
       by the length (int32_t) of the array, followed by the values.
    */
    auto it = std::ranges::begin(stream_view);
    uint16_t tag = static_cast<uint16_t>(*it) << 8;
    ++it; // skip char read before

    tag += static_cast<uint16_t>(*it);
    ++it; // skip char read before

    char type_id = *it;
    ++it; // skip char read before

    switch (type_id)
    {
        case 'A' : // char
        {
            target[tag] = *it;
            ++it; // skip char that has been read
            break;
        }
        // all integer sizes are possible
        case 'c' : // int8_t
        {
            int8_t tmp;
            read_field(stream_view, tmp);
            target[tag] = static_cast<int32_t>(tmp); // readable sam format only allows int32_t
            break;
        }
        case 'C' : // uint8_t
        {
            uint8_t tmp;
            read_field(stream_view, tmp);
            target[tag] = static_cast<int32_t>(tmp); // readable sam format only allows int32_t
            break;
        }
        case 's' : // int16_t
        {
            int16_t tmp;
            read_field(stream_view, tmp);
            target[tag] = static_cast<int32_t>(tmp); // readable sam format only allows int32_t
            break;
        }
        case 'S' : // uint16_t
        {
            uint16_t tmp;
            read_field(stream_view, tmp);
            target[tag] = static_cast<int32_t>(tmp); // readable sam format only allows int32_t
            break;
        }
        case 'i' : // int32_t
        {
            int32_t tmp;
            read_field(stream_view, tmp);
            target[tag] = std::move(tmp); // readable sam format only allows int32_t
            break;
        }
        case 'I' : // uint32_t
        {
            uint32_t tmp;
            read_field(stream_view, tmp);
            target[tag] = static_cast<int32_t>(tmp); // readable sam format only allows int32_t
            break;
        }
        case 'f' : // float
        {
            float tmp;
            read_field(stream_view, tmp);
            target[tag] = tmp;
            break;
        }
        case 'Z' : // string
        {
            string_buffer.clear();
            while (!is_char<'\0'>(*it))
            {
                string_buffer.push_back(*it);
                ++it;
            }
            ++it; // skip \0
            target[tag] = string_buffer;
            break;
        }
        case 'H' : // byte array, represented as null-terminated string; specification requires even number of bytes
        {
            std::vector<std::byte> byte_array;
            std::byte value;
            while (!is_char<'\0'>(*it))
            {
                string_buffer.clear();
                string_buffer.push_back(*it);
                ++it;

                if (*it == '\0')
                    throw format_error{"Hexadecimal tag has an uneven number of digits!"};

                string_buffer.push_back(*it);
                ++it;
                read_field(string_buffer, value);
                byte_array.push_back(value);
            }
            ++it; // skip \0
            target[tag] = byte_array;
            break;
        }
        case 'B' : // Array. Value type depends on second char [cCsSiIf]
        {
            char array_value_type_id = *it;
            ++it; // skip char read before

            switch (array_value_type_id)
            {
                case 'c' : // int8_t
                    read_sam_dict_vector(target[tag], stream_view, int8_t{});
                    break;
                case 'C' : // uint8_t
                    read_sam_dict_vector(target[tag], stream_view, uint8_t{});
                    break;
                case 's' : // int16_t
                    read_sam_dict_vector(target[tag], stream_view, int16_t{});
                    break;
                case 'S' : // uint16_t
                    read_sam_dict_vector(target[tag], stream_view, uint16_t{});
                    break;
                case 'i' : // int32_t
                    read_sam_dict_vector(target[tag], stream_view, int32_t{});
                    break;
                case 'I' : // uint32_t
                    read_sam_dict_vector(target[tag], stream_view, uint32_t{});
                    break;
                case 'f' : // float
                    read_sam_dict_vector(target[tag], stream_view, float{});
                    break;
                default:
                    throw format_error{detail::to_string("The first character in the numerical id of a SAM tag ",
                                       "must be one of [cCsSiIf] but '", array_value_type_id, "' was given.")};
            }
            break;
        }
        default:
            throw format_error{detail::to_string("The second character in the numerical id of a "
                               "SAM tag must be one of [A,i,Z,H,B,f] but '", type_id, "' was given.")};
    }
}

/*!\brief Parses a cigar string into a vector of operation-count pairs (e.g. (M, 3)).
 * \tparam cigar_input_type The type of a single pass input view over the cigar string; must model
 *                          std::ranges::input_range.
 * \param[in] cigar_input The single pass input view over the cigar string to parse.
 * \param[in] n_cigar_op  The number of cigar elements to read from the cigar_input.
 *
 * \returns A tuple of size three containing (1) std::vector over seqan3::cigar, that describes
 *          the alignment, (2) the aligned reference length, (3) the aligned query sequence length.
 *
 * \details
 *
 * For example, the view over the cigar string "1H4M1D2M2S" will return
 * `{[(H,1), (M,4), (D,1), (M,2), (S,2)], 7, 6}`.
 */
template <typename cigar_input_type>
inline auto format_bam::parse_binary_cigar(cigar_input_type && cigar_input, uint16_t n_cigar_op) const
{
    std::vector<cigar> operations{};
    char operation{'\0'};
    uint32_t count{};
    int32_t ref_length{}, seq_length{};
    uint32_t operation_and_count{}; // in BAM operation and count values are stored within one 32 bit integer
    constexpr char const * cigar_mapping = "MIDNSHP=X*******";
    constexpr uint32_t cigar_mask = 0x0f; // 0000000000001111

    if (n_cigar_op == 0) // [[unlikely]]
        return std::tuple{operations, ref_length, seq_length};

    // parse the rest of the cigar
    // -------------------------------------------------------------------------------------------------------------
    while (n_cigar_op > 0) // until stream is not empty
    {
        std::ranges::copy_n(std::ranges::begin(cigar_input),
                            sizeof(operation_and_count),
                            reinterpret_cast<char*>(&operation_and_count));
        operation = cigar_mapping[operation_and_count & cigar_mask];
        count = operation_and_count >> 4;

        detail::update_alignment_lengths(ref_length, seq_length, operation, count);
        operations.emplace_back(count, cigar::operation{}.assign_char(operation));
        --n_cigar_op;
    }

    return std::tuple{operations, ref_length, seq_length};
}

/*!\brief Writes the optional fields of the seqan3::sam_tag_dictionary.
 * \param[in] tag_dict The tag dictionary to print.
 */
inline std::string format_bam::get_tag_dict_str(sam_tag_dictionary const & tag_dict)
{
    std::string result{};

    auto stream_variant_fn = [&result] (auto && arg) // helper to print an std::variant
    {
        // T is either char, int32_t, float, std::string, or a std::vector<some int>
        using T = std::remove_cvref_t<decltype(arg)>;

        if constexpr (std::same_as<T, int32_t>)
        {
            // always choose the smallest possible representation [cCsSiI]
            size_t const absolute_arg = std::abs(arg);
            auto n = std::countr_zero(std::bit_ceil(absolute_arg + 1u) >> 1u) / 8u;
            bool const negative = arg < 0;
            n = n * n + 2 * negative; // for switch case order

            switch (n)
            {
                case 0:
                {
                    result[result.size() - 1] = 'C';
                    result.append(reinterpret_cast<char const *>(&arg), 1);
                    break;
                }
                case 1:
                {
                    result[result.size() - 1] = 'S';
                    result.append(reinterpret_cast<char const *>(&arg), 2);
                    break;
                }
                case 2:
                {
                    result[result.size() - 1] = 'c';
                    int8_t tmp = static_cast<int8_t>(arg);
                    result.append(reinterpret_cast<char const *>(&tmp), 1);
                    break;
                }
                case 3:
                {
                    result[result.size() - 1] = 's';
                    int16_t tmp = static_cast<int16_t>(arg);
                    result.append(reinterpret_cast<char const *>(&tmp), 2);
                    break;
                }
                default:
                {
                    result.append(reinterpret_cast<char const *>(&arg), 4); // always i
                    break;
                }
            }
        }
        else if constexpr (std::same_as<T, std::string>)
        {
            result.append(reinterpret_cast<char const *>(arg.data()), arg.size() + 1/*+ null character*/);
        }
        else if constexpr (!std::ranges::range<T>) // char, float
        {
            result.append(reinterpret_cast<char const *>(&arg), sizeof(arg));
        }
        else // std::vector of some arithmetic_type type
        {
            int32_t sz{static_cast<int32_t>(arg.size())};
            result.append(reinterpret_cast<char *>(&sz), 4);
            result.append(reinterpret_cast<char const *>(arg.data()),
                          arg.size() * sizeof(std::ranges::range_value_t<T>));
        }
    };

    for (auto & [tag, variant] : tag_dict)
    {
        result.push_back(static_cast<char>(tag / 256));
        result.push_back(static_cast<char>(tag % 256));

        result.push_back(detail::sam_tag_type_char[variant.index()]);

        if (!is_char<'\0'>(detail::sam_tag_type_char_extra[variant.index()]))
            result.push_back(detail::sam_tag_type_char_extra[variant.index()]);

        std::visit(stream_variant_fn, variant);
    }

    return result;
}

} // namespace seqan3
