// -----------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2019, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2019, Knut Reinert & MPI für molekulare Genetik
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at: https://github.com/seqan/seqan3/blob/master/LICENSE.md
// -----------------------------------------------------------------------------------------------------

/*!\file
 * \brief Provides the seqan3::format_sam tag and the seqan3::alignment_file_input_format and
 *        seqan3::alignment_file_output_format specialisation for this tag.
 * \author Svenja Mehringer <svenja.mehringer AT fu-berlin.de>
 */

#pragma once

#include <iterator>
#include <string>
#include <vector>

#include <seqan3/core/char_operations/predicate.hpp>
#include <seqan3/core/concept/core_language.hpp>
#include <seqan3/core/concept/tuple.hpp>
#include <seqan3/core/detail/reflection.hpp>
#include <seqan3/core/detail/to_string.hpp>
#include <seqan3/core/type_traits/range.hpp>
#include <seqan3/core/type_traits/template_inspection.hpp>
#include <seqan3/io/alignment_file/detail.hpp>
#include <seqan3/io/alignment_file/header.hpp>
#include <seqan3/io/alignment_file/input_format_concept.hpp>
#include <seqan3/io/alignment_file/input_options.hpp>
#include <seqan3/io/alignment_file/output_format_concept.hpp>
#include <seqan3/io/alignment_file/output_options.hpp>
#include <seqan3/io/alignment_file/sam_tag_dictionary.hpp>
#include <seqan3/io/detail/ignore_output_iterator.hpp>
#include <seqan3/io/detail/misc.hpp>
#include <seqan3/io/stream/iterator.hpp>
#include <seqan3/range/decorator/gap_decorator.hpp>
#include <seqan3/range/detail/misc.hpp>
#include <seqan3/range/view/char_to.hpp>
#include <seqan3/range/view/istreambuf.hpp>
#include <seqan3/range/view/repeat_n.hpp>
#include <seqan3/range/view/slice.hpp>
#include <seqan3/range/view/take_until.hpp>
#include <seqan3/range/view/to_char.hpp>
#include <seqan3/std/algorithm>
#include <seqan3/std/charconv>
#include <seqan3/std/concepts>
#include <seqan3/std/ranges>
#include <seqan3/version.hpp>

namespace seqan3
{

/*!\brief       The SAM format (tag).
 * \implements  AlignmentFileFormat
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
 * Take a look at our tutorial \ref tutorial_alignment_file for a walk through of how to read alignment files.
 *
 * ### Fields
 *
 * The SAM format provides the following fields:
 * seqan3::field::ALIGNMENT, seqan3::field::SEQ, seqan3::field::QUAL,
 * seqan3::field::ID, seqan3::field::REF_SEQ, seqan3::field::REF_ID
 * seqan3::field::REF_OSSFET, seqan3::field::OFFSET, seqan3::field::FLAG,
 * seqan3::field::MAPQ and seqan3::field::MATE.
 * In addition there is the seqan3::field::HEADER_PTR, which is usually only used internally
 * to provide the range-based functionality of the file.
 *
 * **None of the fields are required** when writing but will be defaulted
 * to '0' for numeric fields and '*' for other fields.
 *
 * ### SAM format columns -> fields
 *
 * Since many users will be accustomed to the columns of the SAM format, here is a
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
 * ### Header implementation
 *
 * The SAM header (if present) is read/written once in the beginning before the
 * first record is read/written.
 */
struct format_sam
{
    //!\brief The valid file extensions for this format; note that you can modify this value.
    static inline std::vector<std::string> file_extensions
    {
        { "sam" },
    };

    //!\brief The format version string.
    static constexpr char format_version[4] = "1.6";
};

} // namespace seqan3

namespace seqan3::detail
{

//!\brief The seqan3::alignment_file_input_format specialisation that handles formatted SAM input.
//!\ingroup alignment_file
template <>
class alignment_file_input_format<format_sam>
{
public:
    //!\brief Exposes the format tag that this class is specialised with
    using format_tag = format_sam;

    /*!\name Constructors, destructor and assignment
     * \{
     */
    alignment_file_input_format()                                                noexcept = default; //!< Defaulted.
    //!\brief Copy construction is explicitly deleted, because you can't have multiple access to the same file.
    alignment_file_input_format(alignment_file_input_format const &)                      = delete;
    //!\brief Copy assignment is explicitly deleted, because you can't have multiple access to the same file.
    alignment_file_input_format & operator=(alignment_file_input_format const &)          = delete;
    alignment_file_input_format(alignment_file_input_format &&)                  noexcept = default; //!< Defaulted.
    alignment_file_input_format & operator=(alignment_file_input_format &&)      noexcept = default; //!< Defaulted.
    ~alignment_file_input_format()                                               noexcept = default; //!< Defaulted.
    //!\}

    //!\copydoc AlignmentFileInputFormat::read
    template <typename stream_type,     // constraints checked by file
              typename seq_legal_alph_type,
              typename ref_seqs_type,
              typename ref_ids_type,
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
    void read(stream_type                                             & stream,
              alignment_file_input_options<seq_legal_alph_type> const & SEQAN3_DOXYGEN_ONLY(options),
              ref_seqs_type                                           & ref_seqs,
              alignment_file_header<ref_ids_type>                     & header,
              seq_type                                                & seq,
              qual_type                                               & qual,
              id_type                                                 & id,
              offset_type                                             & offset,
              ref_seq_type                                            & SEQAN3_DOXYGEN_ONLY(ref_seq),
              ref_id_type                                             & ref_id,
              ref_offset_type                                         & ref_offset,
              align_type                                              & align,
              flag_type                                               & flag,
              mapq_type                                               & mapq,
              mate_type                                               & mate,
              tag_dict_type                                           & tag_dict,
              e_value_type                                            & SEQAN3_DOXYGEN_ONLY(e_value),
              bit_score_type                                          & SEQAN3_DOXYGEN_ONLY(bit_score))
    {
        static_assert(detail::decays_to_ignore_v<ref_offset_type> ||
                      detail::is_type_specialisation_of_v<ref_offset_type, std::optional>,
                      "The ref_offset must be a specialisation of std::optional.");

        auto stream_view = view::istreambuf(stream);
        auto field_view = stream_view | view::take_until_or_throw_and_consume(is_char<'\t'>);

        // these variables need to be stored to compute the ALIGNMENT
        int32_t ref_offset_tmp{};
        value_type_t<decltype(header.ref_ids())> ref_id_tmp{};
        [[maybe_unused]] int32_t offset_tmp{};
        [[maybe_unused]] int32_t soft_clipping_end{};
        [[maybe_unused]] std::vector<std::pair<char, size_t>> cigar{};
        [[maybe_unused]] int32_t ref_length{0}, seq_length{0}; // length of aligned part for ref and query

        // Header
        // -------------------------------------------------------------------------------------------------------------
        if (is_char<'@'>(*std::ranges::begin(stream_view))) // we always read the header if present
        {
            read_header(stream_view, header, ref_seqs);

            if (std::ranges::begin(stream_view) == std::ranges::end(stream_view)) // file has no records
                return;
        }

        // Fields 1-5: ID FLAG REF_ID REF_OFFSET MAPQ
        // -------------------------------------------------------------------------------------------------------------
        read_field(field_view, id);

        read_field(field_view, flag);

        read_field(field_view, ref_id_tmp);
        check_and_assign_ref_id(ref_id, ref_id_tmp, header, ref_seqs);

        read_field(field_view, ref_offset_tmp);
        --ref_offset_tmp; // SAM format is 1-based but SeqAn operates 0-based

        if (ref_offset_tmp == -1)
            ref_offset = std::nullopt; // indicates an unmapped read -> ref_offset is not set
        else if (ref_offset_tmp > -1)
            ref_offset = ref_offset_tmp;
        else if (ref_offset_tmp < -1)
            throw format_error{"No negative values are allowed for field::REF_OFFSET."};

        read_field(field_view, mapq);

        // Field 6: CIGAR
        // -------------------------------------------------------------------------------------------------------------
        if constexpr (!detail::decays_to_ignore_v<align_type>)
        {
            if (!is_char<'*'>(*std::ranges::begin(stream_view))) // no cigar information given
            {
                std::tie(cigar, ref_length, seq_length, offset_tmp, soft_clipping_end) = parse_cigar(field_view);
            }
            else
            {
                std::ranges::next(std::ranges::begin(field_view)); // skip '*'
            }
        }
        else
        {
            detail::consume(field_view);
        }

        offset = offset_tmp;

        // Field 7-9: (RNEXT PNEXT TLEN) = MATE
        // -------------------------------------------------------------------------------------------------------------
        if constexpr (!detail::decays_to_ignore_v<mate_type>)
        {
            value_type_t<decltype(header.ref_ids())> tmp_mate_ref_id{};
            read_field(field_view, tmp_mate_ref_id); // RNEXT

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
            read_field(field_view, tmp_pnext); // PNEXT

            if (tmp_pnext > 0)
                get<1>(mate) = --tmp_pnext; // SAM format is 1-based but SeqAn operates 0-based.
            else if (tmp_pnext < 0)
                throw format_error{"No negative values are allowed at the mate mapping position."};
            // tmp_pnext == 0 indicates an unmapped mate -> do not fill std::optional get<1>(mate)

            read_field(field_view, get<2>(mate)); // TLEN
        }
        else
        {
            for (size_t i = 0; i < 3u; ++i)
            {
                detail::consume(field_view);
            }
        }

        // Field 10: Sequence
        // -------------------------------------------------------------------------------------------------------------
        if (!is_char<'*'>(*std::ranges::begin(stream_view))) // sequence information is given
        {
            auto constexpr is_legal_alph = is_in_alphabet<seq_legal_alph_type>;
            auto seq_stream = field_view | std::view::transform([is_legal_alph] (char const c) // enforce legal alphabet
                                           {
                                               if (!is_legal_alph(c))
                                                   throw format_error{std::string{"Encountered an unexpected letter: "} +
                                                                     is_legal_alph.msg.str() +
                                                                     " evaluated to false on " +
                                                                     detail::make_printable(c)};
                                               return c;
                                           });

            if constexpr (detail::decays_to_ignore_v<seq_type>)
            {
                if constexpr (!detail::decays_to_ignore_v<align_type>)
                {
                    static_assert(SequenceContainer<std::remove_reference_t<decltype(get<1>(align))>>,
                                  "If you want to read ALIGNMENT but not SEQ, the alignment"
                                  " object must store a sequence container at the second (query) position.");

                    if (!cigar.empty()) // only parse alignment if cigar information was given
                    {

                        auto tmp_iter = std::ranges::begin(seq_stream);
                        std::ranges::advance(tmp_iter, offset_tmp);

                        for (; seq_length > 0; --seq_length) // seq_length is not needed anymore
                        {
                            get<1>(align).push_back(value_type_t<decltype(get<1>(align))>{}.assign_char(*tmp_iter));
                            ++tmp_iter;
                        }

                        std::ranges::advance(tmp_iter, soft_clipping_end);
                    }
                    else
                    {
                        get<1>(align) = std::remove_reference_t<decltype(get<1>(align))>{}; // empty container
                    }
                }
                else
                {
                    detail::consume(seq_stream);
                }
            }
            else
            {
                read_field(seq_stream, seq);

                if constexpr (!detail::decays_to_ignore_v<align_type>)
                {
                    if (!cigar.empty()) // if no alignment info is given, the field::ALIGNMENT should remain empty
                    {
                        assign_unaligned(get<1>(align),
                                         seq | view::slice(static_cast<decltype(std::ranges::size(seq))>(offset_tmp),
                                                           std::ranges::size(seq) - soft_clipping_end));
                    }
                }
            }
        }
        else
        {
            std::ranges::next(std::ranges::begin(field_view)); // skip '*'
        }

        // Field 11:  Quality
        // -------------------------------------------------------------------------------------------------------------
        auto const tab_or_end = is_char<'\t'> || is_char<'\r'> || is_char<'\n'>;
        read_field(stream_view | view::take_until_or_throw(tab_or_end), qual);

        if constexpr (!detail::decays_to_ignore_v<seq_type> && !detail::decays_to_ignore_v<qual_type>)
        {
            if (std::ranges::distance(seq) != 0 && std::ranges::distance(qual) != 0 &&
                std::ranges::distance(seq) != std::ranges::distance(qual))
            {
                throw format_error{to_string("Sequence length (", std::ranges::distance(seq), ") and quality length (",
                                             std::ranges::distance(qual), ") must be the same.")};
            }
        }

        // All remaining optional fields if any: SAM tags dictionary
        // -------------------------------------------------------------------------------------------------------------
        while (is_char<'\t'>(*std::ranges::begin(stream_view))) // read all tags if present
        {
            std::ranges::next(std::ranges::begin(stream_view)); // skip tab
            read_field(stream_view | view::take_until_or_throw(tab_or_end), tag_dict);
        }

        detail::consume(stream_view | view::take_until(!(is_char<'\r'> || is_char<'\n'>))); // consume new line

        // DONE READING - wrap up
        // -------------------------------------------------------------------------------------------------------------
        // Alignment object construction
        // Note that the query sequence in get<1>(align) has already been filled while reading Field 10.
        if constexpr (!detail::decays_to_ignore_v<align_type>)
        {
            int32_t ref_idx{(ref_id_tmp.empty()/*unmapped read?*/) ? -1 : 0};

            if constexpr (!detail::decays_to_ignore_v<ref_seqs_type>)
            {
                if (!ref_id_tmp.empty())
                {
                    assert(header.ref_dict.count(ref_id_tmp) != 0); // taken care of in check_and_assign_ref_id()
                    ref_idx = header.ref_dict[ref_id_tmp];          // get index for reference sequence
                }
            }

            construct_alignment(align, cigar, ref_idx, ref_seqs, ref_offset_tmp, ref_length);
        }
    }

protected:
    //!\privatesection
    //!\brief A buffer used when parsing arithmetic values with std::from_chars.
    std::array<char, 316> buffer{}; // Doubles can be up to 316 characters

    //!\brief Tracks whether reference information (\@SR tag) were found in the header
    bool ref_info_present_in_header{false};

    /*!\brief Checks for known reference ids or adds a new reference is and assigns a reference id to `ref_id`.
     * \tparam ref_id_type         The type of the reference id (usually a view::all over ref_id_tmp_type).
     * \tparam ref_id_tmp_type     The type of the temporary parsed id (Same type as reference ids in header).
     * \tparam header_type         The type of the alignment header.
     * \tparam ref_seqs_type       A tag whether the reference information were given or not (std::ignore or not).
     *
     * \param[in, out] ref_id      The reference id to be filled.
     * \param[in, out] ref_id_tmp  The temporary of the parsed reference id.
     * \param[in, out] header      The header object that stores the reference id information.
     */
    template <typename ref_id_type,
              typename ref_id_tmp_type,
              typename header_type,
              typename ref_seqs_type>
    void check_and_assign_ref_id(ref_id_type      & ref_id,
                                 ref_id_tmp_type  & ref_id_tmp,
                                 header_type      & header,
                                 ref_seqs_type    & /*tag*/)
    {
        if (!std::ranges::empty(ref_id_tmp)) // otherwise the std::optional will not be filled
        {
            auto search = header.ref_dict.find(ref_id_tmp);

            if (search == header.ref_dict.end())
            {
                if constexpr(detail::decays_to_ignore_v<ref_seqs_type>) // no reference information given
                {
                    if (ref_info_present_in_header)
                    {
                        throw format_error{"Unknown reference id found in record which is not present in the header."};
                    }
                    else
                    {
                        header.ref_ids().push_back(ref_id_tmp);
                        auto pos = std::ranges::size(header.ref_ids()) - 1;
                        header.ref_dict[header.ref_ids()[pos]] = pos;
                        ref_id = pos;
                    }
                }
                else
                {
                    throw format_error{"Unknown reference id found in record which is not present in the given ids."};
                }
            }
            else
            {
                ref_id = search->second;
            }
        }
    }

    /*!\brief Construct the field::ALIGNMENT depending on the given information.
     * \tparam align_type    The alignment type.
     * \tparam ref_seqs_type The type of reference sequences (might decay to ignore).
     * \param[in,out] align  The alignment (pair of aligned sequences) to fill.
     * \param[in] cigar      The cigar information to convert to an alignment.
     * \param[in] rid        The index of the reference sequence in header.ref_ids().
     * \param[in] ref_seqs   The reference sequence information.
     * \param[in] ref_start  The start position of the alignment in the reference sequence.
     * \param[in] ref_length The length of the aligned reference sequence.
     */
    template <typename align_type, typename ref_seqs_type>
    void construct_alignment(align_type                           & align,
                             std::vector<std::pair<char, size_t>> & cigar,
                             [[maybe_unused]] int32_t               rid,
                             [[maybe_unused]] ref_seqs_type       & ref_seqs,
                             [[maybe_unused]] int32_t               ref_start,
                             size_t                                 ref_length)
    {
        if (rid > -1 && ref_start > -1 &&       // read is mapped
            !cigar.empty() &&                   // alignment field was not empty
            !std::ranges::empty(get<1>(align))) // seq field was not empty
        {
            if constexpr (!detail::decays_to_ignore_v<ref_seqs_type>)
            {
                assert(static_cast<size_t>(ref_start + ref_length) <= std::ranges::size(ref_seqs[rid]));
                // copy over unaligned reference sequence part
                assign_unaligned(get<0>(align), ref_seqs[rid] | view::slice(ref_start, ref_start + ref_length));
            }
            else
            {
                using unaligned_t = remove_cvref_t<detail::unaligned_seq_t<decltype(get<0>(align))>>;
                auto dummy_seq    = view::repeat_n(value_type_t<unaligned_t>{}, ref_length)
                                  | std::view::transform(detail::access_restrictor_fn{});
                static_assert(std::Same<unaligned_t, decltype(dummy_seq)>,
                              "No reference information was given so the type of the first alignment tuple position"
                              "must have an unaligned sequence type of a dummy sequence ("
                              "view::repeat_n(dna5{}, size_t{}) | "
                              "std::view::transform(detail::access_restrictor_fn{}))");

                assign_unaligned(get<0>(align), dummy_seq); // assign dummy sequence
            }

            // insert gaps according to the cigar information
            detail::alignment_from_cigar(align, cigar);
        }
        else // not enough information for an alignment, assign an empty view/dummy_sequence
        {
            if constexpr (!detail::decays_to_ignore_v<ref_seqs_type>) // reference info given
            {
                assert(std::ranges::size(ref_seqs) > 0); // we assume that the given ref info is not empty
                assign_unaligned(get<0>(align), ref_seqs[0] | view::slice(0, 0));
            }
            else
            {
                using unaligned_t = remove_cvref_t<detail::unaligned_seq_t<decltype(get<0>(align))>>;
                assign_unaligned(get<0>(align), view::repeat_n(value_type_t<unaligned_t>{}, 0)
                                                | std::view::transform(detail::access_restrictor_fn{}));
            }
        }
    }

    /*!\brief Decays to detail::consume for std::ignore.
     * \tparam stream_view_type  The type of the stream as a view.
     *
     * \param[in, out] stream_view  The stream view to consume.
     * \param[in, out] target       A std::ignore placeholder.
     */
    template <typename stream_view_type>
    void read_field(stream_view_type && stream_view, detail::ignore_t const & SEQAN3_DOXYGEN_ONLY(target))
    {
        detail::consume(stream_view);
    }

    /*!\brief Reads a range by copying from stream_view to target, converting values with seqan3::view::char_to.
     * \tparam stream_view_type  The type of the stream as a view.
     * \tparam target_range_type The type of range to parse from input; must model std::ranges::ForwardRange.
     *
     * \param[in, out] stream_view  The stream view to iterate over.
     * \param[in, out] target       The range to store the parsed sequence.
     */
    template <typename stream_view_type, std::ranges::ForwardRange target_range_type>
    void read_field(stream_view_type && stream_view, target_range_type & target)
    {
        if (!is_char<'*'>(*std::ranges::begin(stream_view)))
            std::ranges::copy(stream_view | view::char_to<value_type_t<target_range_type>>,
                              std::back_inserter(target));
        else
            std::ranges::next(std::ranges::begin(stream_view)); // skip '*'
    }

    /*!\brief Reads arithmetic fields using std::from_chars.
     * \tparam stream_view_type The type of the stream as a view.
     * \tparam target_type      The type of value to parse from input; must model seqan3::Arithmetic.
     *
     * \param[in, out] stream_view  The stream view to iterate over.
     * \param[in, out] target       The arithmetic value object to store the parsed value.
     *
     * \throws seqan3::format_error if the character sequence in stream_view cannot be successfully converted to a value
     *         of type target_type.
     */
    template <typename stream_view_type, Arithmetic target_type>
    void read_field(stream_view_type && stream_view, target_type & target)
    {
        // unfortunately std::from_chars only accepts char const * so we need a buffer.
        auto [ignore, end] = std::ranges::copy(stream_view, buffer.data());
        (void) ignore;
        std::from_chars_result res = std::from_chars(buffer.begin(), end, target);

        if (res.ec == std::errc::invalid_argument || res.ptr != end)
            throw format_error{std::string("[CORRUPTED SAM FILE] The string '") + std::string(buffer.begin(), end) +
                                           "' could not be cast into type " +
                                           detail::get_display_name_v<target_type>.str()};

        if (res.ec == std::errc::result_out_of_range)
            throw format_error{std::string("[CORRUPTED SAM FILE] Casting '") + std::string(buffer.begin(), end) +
                                           "' into type " + detail::get_display_name_v<target_type>.str() +
                                           " would cause an overflow."};
    }

    /*!\brief Delegate parsing of std::optional types to parsing of the inner value type.
     * \tparam stream_view_type     The type of the stream as a view.
     * \tparam optional_value_type  The inner type of a the std::optional type of \p target.
     *
     * \param[in, out] stream_view  The stream view to iterate over.
     * \param[in, out] target       The std::optional object to store the parsed value.
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

    /*!\brief Reads a list of values separated by comma as it is the case for SAM tag arrays.
     * \tparam stream_view_type The type of the stream as a view.
     * \tparam value_type       The type of values to be stored in the tag array.
     *
     * \param[in, out] variant      A std::variant object to store the tag arrays.
     * \param[in, out] stream_view  The stream view to iterate over.
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
    template <typename stream_view_type, typename value_type>
    void read_sam_dict_vector(seqan3::detail::sam_tag_variant & variant,
                              stream_view_type && stream_view,
                              value_type value)
    {
        std::vector<value_type> tmp_vector;
        while (std::ranges::begin(stream_view) != ranges::end(stream_view)) // not fully consumed yet
        {
            read_field(stream_view | view::take_until(is_char<','>), value);
            tmp_vector.push_back(value);

            if (is_char<','>(*std::ranges::begin(stream_view)))
                std::ranges::next(std::ranges::begin(stream_view)); // skip ','
        }
        variant = std::move(tmp_vector);
    }

    /*!\brief Reads the optional tag fields into the seqan3::sam_tag_dictionary.
     * \tparam stream_view_type   The type of the stream as a view.
     *
     * \param[in, out] stream_view  The stream view to iterate over.
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
    template <typename stream_view_type>
    void read_field(stream_view_type && stream_view, sam_tag_dictionary & target)
    {
        /* Every SAM tag has the format "[TAG]:[TYPE_ID]:[VALUE]", where TAG is a two letter
           name tag which is converted to a unique integer identifier and TYPE_ID is one character in [A,i,Z,H,B,f]
           describing the type for the upcoming VALUES. If TYPE_ID=='B' it signals an array of comma separated
           VALUE's and the inner value type is identified by the character following ':', one of [cCsSiIf].
        */
        uint16_t tag = static_cast<uint16_t>(*std::ranges::begin(stream_view)) << 8;
        std::ranges::next(std::ranges::begin(stream_view)); // skip char read before
        tag += static_cast<uint16_t>(*std::ranges::begin(stream_view));
        std::ranges::next(std::ranges::begin(stream_view)); // skip char read before
        std::ranges::next(std::ranges::begin(stream_view)); // skip ':'
        char type_id = *std::ranges::begin(stream_view);
        std::ranges::next(std::ranges::begin(stream_view)); // skip char read before
        std::ranges::next(std::ranges::begin(stream_view)); // skip ':'

        switch (type_id)
        {
            case 'A' : // char
            {
                target[tag] = static_cast<char>(*std::ranges::begin(stream_view));
                std::ranges::next(std::ranges::begin(stream_view)); // skip char that has been read
                break;
            }
            case 'i' : // int32_t
            {
                int32_t tmp;
                read_field(stream_view, tmp);
                target[tag] = tmp;
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
                target[tag] = std::string(stream_view);
                break;
            }
            case 'H' :
            {
                // TODO
                break;
            }
            case 'B' : // Array. Value type depends on second char [cCsSiIf]
            {
                char array_value_type_id = *std::ranges::begin(stream_view);
                std::ranges::next(std::ranges::begin(stream_view)); // skip char read before
                std::ranges::next(std::ranges::begin(stream_view)); // skip first ','

                switch (array_value_type_id)
                {
                    case 'c' : // int8_t
                    {
                        read_sam_dict_vector(target[tag], stream_view, int8_t{});
                        break;
                    }
                    case 'C' : // uint8_t
                    {
                        read_sam_dict_vector(target[tag], stream_view, uint8_t{});
                        break;
                    }
                    case 's' : // int16_t
                    {
                        read_sam_dict_vector(target[tag], stream_view, int16_t{});
                        break;
                    }
                    case 'S' : // uint16_t
                    {
                        read_sam_dict_vector(target[tag], stream_view, uint16_t{});
                        break;
                    }
                    case 'i' : // int32_t
                    {
                        read_sam_dict_vector(target[tag], stream_view, int32_t{});
                        break;
                    }
                    case 'I' : // uint32_t
                    {
                        read_sam_dict_vector(target[tag], stream_view, uint32_t{});
                        break;
                    }
                    case 'f' : // float
                    {
                        read_sam_dict_vector(target[tag], stream_view, float{});
                        break;
                    }
                    default:
                        throw format_error{std::string("The first character in the numerical ") +
                                           "id of a SAM tag must be one of [cCsSiIf] but '" + array_value_type_id +
                                           "' was given."};
                }
                break;
            }
            default:
                throw format_error{std::string("The second character in the numerical id of a "
                                   "SAM tag must be one of [A,i,Z,H,B,f] but '") + type_id + "' was given."};
        }
    }

    /*!\brief Reads the SAM header.
     * \tparam stream_view_type     The type of the stream as a view.
     * \param[in, out] stream_view  The stream view to iterate over.
     * \param[in, out] hdr          The header (as a pointer) to store the parsed values.
     *
     * \throws seqan3::format_error if any unexpected character or format is encountered.
     *
     * \details
     *
     * Reading the header format is done according to the official
     * [SAM format specifications](https://samtools.github.io/hts-specs/SAMv1.pdf).
     *
     * The function throws a seqan3::format_error if any unknown tag was encountered. It will also fail if the format is
     * not in a correct state (e.g. required fields are not given), but throwing might occur downstream of the actual
     * error.
     */
    template <typename stream_view_type, typename ref_ids_type, typename ref_seqs_type>
    void read_header(stream_view_type && stream_view,
                     alignment_file_header<ref_ids_type> & hdr,
                     ref_seqs_type & /*ref_id_to_pos_map*/)
    {
        auto parse_tag_value = [&stream_view, this] (auto & value) // helper function to parse the next tag value
        {
            detail::consume(stream_view | view::take_until_or_throw(is_char<':'>)); // skip tag name
            std::ranges::next(std::ranges::begin(stream_view));                     // skip ':'
            read_field(stream_view | view::take_until_or_throw(is_char<'\t'> || is_char<'\n'>), value);
        };

        // @HQ line
        // -------------------------------------------------------------------------------------------------------------
        parse_tag_value(hdr.format_version); // parse required VN (version) tag

        // The SO, SS and GO tag are optional and can appear in any order
        while (is_char<'\t'>(*std::ranges::begin(stream_view)))
        {
            std::ranges::next(std::ranges::begin(stream_view));              // skip tab
            std::string * who = std::addressof(hdr.grouping);

            if (is_char<'S'>(*std::ranges::begin(stream_view)))
            {
                std::ranges::next(std::ranges::begin(stream_view));          // skip S

                if (is_char<'O'>(*std::ranges::begin(stream_view)))          // SO (sorting) tag
                    who = std::addressof(hdr.sorting);
                else if (is_char<'S'>(*std::ranges::begin(stream_view)))     // SS (sub-order) tag
                    who = std::addressof(hdr.subsorting);
                else
                    throw format_error{std::string{"Illegal SAM header tag: S"} +
                                       std::string{static_cast<char>(*std::ranges::begin(stream_view))}};
            }
            else if (!is_char<'G'>(*std::ranges::begin(stream_view)))        // GO (grouping) tag
            {
                throw format_error{std::string{"Illegal SAM header tag in @HG starting with:"} +
                                   std::string{static_cast<char>(*std::ranges::begin(stream_view))}};
            }

            parse_tag_value(*who);
        }
        std::ranges::next(std::ranges::begin(stream_view));                  // skip newline

        // The rest of the header lines
        // -------------------------------------------------------------------------------------------------------------
        while (is_char<'@'>(*std::ranges::begin(stream_view)))
        {
            std::ranges::next(std::ranges::begin(stream_view));              // skip @

            if (is_char<'S'>(*std::ranges::begin(stream_view)))              // SQ (sequence dictionary) tag
            {
                ref_info_present_in_header = true;
                value_type_t<decltype(hdr.ref_ids())> id;
                std::tuple<int32_t, std::string> info{};

                parse_tag_value(id);                                         // parse required SN (sequence name) tag
                std::ranges::next(std::ranges::begin(stream_view));          // skip tab or newline
                parse_tag_value(get<0>(info));                               // parse required LN (length) tag

                if (is_char<'\t'>(*std::ranges::begin(stream_view)))         // read rest of the tags
                {
                    std::ranges::next(std::ranges::begin(stream_view));      // skip tab
                    read_field(stream_view | view::take_until_or_throw(is_char<'\n'>), get<1>(info));
                }
                std::ranges::next(std::ranges::begin(stream_view));          // skip newline

                /* If reference information were given, the ids exist and we can fill ref_dict directly.
                 * If not, wee need to update the ids first and fill the reference dictionary afterwards. */
                if constexpr (!detail::decays_to_ignore_v<ref_seqs_type>) // reference information given
                {
                    auto id_it = hdr.ref_dict.find(id);

                    if (id_it == hdr.ref_dict.end())
                        throw format_error{to_string("Unknown reference name '", id, "' found in SAM header ",
                                                     "(header.ref_ids(): ", hdr.ref_ids(), ").")};

                    auto & given_ref_info = hdr.ref_id_info[id_it->second];

                    if (std::get<0>(given_ref_info) != std::get<0>(info))
                        throw format_error{"Provided reference has unequal length as specified in the header."};

                    hdr.ref_id_info[id_it->second] = std::move(info);
                }
                else
                {
                    static_assert(!detail::is_type_specialisation_of_v<decltype(hdr.ref_ids()), std::deque>,
                                  "The range over reference ids must be of type std::deque such that "
                                  "pointers are not invalidated.");

                    hdr.ref_ids().push_back(id);
                    hdr.ref_id_info.push_back(info);
                    hdr.ref_dict[(hdr.ref_ids())[(hdr.ref_ids()).size() - 1]] = (hdr.ref_ids()).size() - 1;
                }
            }
            else if (is_char<'R'>(*std::ranges::begin(stream_view)))         // RG (read group) tag
            {
                std::pair<std::string, std::string> tmp{};

                parse_tag_value(get<0>(tmp));                                // read required ID tag

                if (is_char<'\t'>(*std::ranges::begin(stream_view)))         // read rest of the tags
                {
                    std::ranges::next(std::ranges::begin(stream_view));
                    read_field(stream_view | view::take_until_or_throw(is_char<'\n'>), get<1>(tmp));
                }
                std::ranges::next(std::ranges::begin(stream_view));          // skip newline

                hdr.read_groups.emplace_back(std::move(tmp));
            }
            else if (is_char<'P'>(*std::ranges::begin(stream_view)))         // PG (program) tag
            {
                typename alignment_file_header<ref_ids_type>::program_info_t tmp{};

                parse_tag_value(tmp.id);                                     // read required ID tag

                // The PN, CL, PP, DS, VN are optional tags and can be given in any order.
                while (is_char<'\t'>(*std::ranges::begin(stream_view)))
                {
                    std::ranges::next(std::ranges::begin(stream_view));      // skip tab
                    std::string * who = &tmp.version;

                    if (is_char<'P'>(*std::ranges::begin(stream_view)))
                    {
                        std::ranges::next(std::ranges::begin(stream_view));  // skip P

                        if (is_char<'N'>(*std::ranges::begin(stream_view)))  // PN (program name) tag
                            who = &tmp.name;
                        else                                                 // PP (previous program) tag
                            who = &tmp.previous;
                    }
                    else if (is_char<'C'>(*std::ranges::begin(stream_view))) // CL (command line) tag
                    {
                        who = &tmp.command_line_call;
                    }
                    else if (is_char<'D'>(*std::ranges::begin(stream_view))) // DS (description) tag
                    {
                        who = &tmp.description;
                    }
                    else if (!is_char<'V'>(*std::ranges::begin(stream_view))) // VN (version) tag
                    {
                        throw format_error{std::string{"Illegal SAM header tag starting with:"} +
                                           std::string{static_cast<char>(*std::ranges::begin(stream_view))}};
                    }

                    parse_tag_value(*who);
                }
                std::ranges::next(std::ranges::begin(stream_view));          // skip newline

                hdr.program_infos.emplace_back(std::move(tmp));
            }
            else if (is_char<'C'>(*std::ranges::begin(stream_view)))         // CO (comment) tag
            {
                std::string tmp;
                std::ranges::next(std::ranges::begin(stream_view)); // skip C
                std::ranges::next(std::ranges::begin(stream_view)); // skip O
                std::ranges::next(std::ranges::begin(stream_view)); // skip :
                read_field(stream_view | view::take_until_or_throw(is_char<'\n'>), tmp);
                std::ranges::next(std::ranges::begin(stream_view)); // skip newline

                hdr.comments.emplace_back(std::move(tmp));
            }
            else
            {
                throw format_error{std::string{"Illegal SAM header tag starting with:"} +
                                   std::string{static_cast<char>(*std::ranges::begin(stream_view))}};
            }
        }
    }
};

//!\brief The seqan3::alignment_file_output_format specialisation that can write formatted SAM.
//!\ingroup alignment_file
template <>
class alignment_file_output_format<format_sam>
{
public:
    //!\brief Exposes the format tag that this class is specialised with
    using format_tag = format_sam;

    /*!\name Constructors, destructor and assignment
     * \{
     */
    alignment_file_output_format()                                                 noexcept = default; //!< Defaulted.
    //!\brief Copy construction is explicitly deleted, because you can't have multiple access to the same file.
    alignment_file_output_format(alignment_file_output_format const &)                      = delete;
    //!\brief Copy assignment is explicitly deleted, because you can't have multiple access to the same file.
    alignment_file_output_format & operator=(alignment_file_output_format const &)          = delete;
    alignment_file_output_format(alignment_file_output_format &&)                  noexcept = default; //!< Defaulted.
    alignment_file_output_format & operator=(alignment_file_output_format &&)      noexcept = default; //!< Defaulted.
    ~alignment_file_output_format()                                                noexcept = default; //!< Defaulted.
    //!\}

    //!\copydoc AlignmentFileOutputFormat::write
    template <typename stream_type,
              typename header_type,
              typename seq_type,
              typename id_type,
              typename ref_seq_type,
              typename ref_id_type,
              typename align_type,
              typename qual_type,
              typename mate_type,
              typename tag_dict_type,
              typename e_value_type,
              typename bit_score_type>
    void write(stream_type                            &  stream,
               alignment_file_output_options const    &  options,
               header_type                            && header,
               seq_type                               && seq,
               qual_type                              && qual,
               id_type                                && id,
               int32_t                                   offset,
               ref_seq_type                           && SEQAN3_DOXYGEN_ONLY(ref_seq),
               ref_id_type                            && ref_id,
               std::optional<int32_t>                    ref_offset,
               align_type                             && align,
               uint16_t                                  flag,
               uint8_t                                   mapq,
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
                      Alphabet<reference_t<seq_type>>),
                      "The seq object must be a std::ranges::ForwardRange over "
                      "letters that model seqan3::Alphabet.");

        static_assert((std::ranges::ForwardRange<id_type>         &&
                      Alphabet<reference_t<id_type>>),
                      "The id object must be a std::ranges::ForwardRange over "
                      "letters that model seqan3::Alphabet.");

        if constexpr (!detail::decays_to_ignore_v<ref_id_type>)
        {
            static_assert((std::ranges::ForwardRange<ref_id_type> ||
                           std::Integral<std::remove_reference_t<ref_id_type>> ||
                           detail::is_type_specialisation_of_v<remove_cvref_t<ref_id_type>, std::optional>),
                          "The ref_id object must be a std::ranges::ForwardRange "
                          "over letters that model seqan3::Alphabet.");

            if constexpr (std::Integral<remove_cvref_t<ref_id_type>> ||
                          detail::is_type_specialisation_of_v<remove_cvref_t<ref_id_type>, std::optional>)
                static_assert(!detail::decays_to_ignore_v<header_type>,
                              "If you give indices as reference id information the header must also be present.");
        }

        static_assert(TupleLike<remove_cvref_t<align_type>>,
                      "The align object must be a std::pair of two ranges whose "
                      "value_type is comparable to seqan3::gap");

        static_assert((std::tuple_size_v<remove_cvref_t<align_type>> == 2 &&
                       std::EqualityComparableWith<gap, reference_t<decltype(std::get<0>(align))>> &&
                       std::EqualityComparableWith<gap, reference_t<decltype(std::get<1>(align))>>),
                      "The align object must be a std::pair of two ranges whose "
                      "value_type is comparable to seqan3::gap");

        static_assert((std::ranges::ForwardRange<qual_type>       &&
                       Alphabet<reference_t<qual_type>>),
                      "The qual object must be a std::ranges::ForwardRange "
                      "over letters that model seqan3::Alphabet.");

        static_assert(TupleLike<remove_cvref_t<mate_type>>,
                      "The mate object must be a std::tuple of size 3 with "
                      "1) a std::ranges::ForwardRange with a value_type modelling seqan3::Alphabet, "
                      "2) a std::Integral or std::optional<std::Integral>, and "
                      "3) a std::Integral.");

        static_assert(((std::ranges::ForwardRange<decltype(std::get<0>(mate))>     ||
                        std::Integral<remove_cvref_t<decltype(std::get<0>(mate))>> ||
                        detail::is_type_specialisation_of_v<remove_cvref_t<decltype(std::get<0>(mate))>, std::optional>) &&
                      (std::Integral<remove_cvref_t<decltype(std::get<1>(mate))>> ||
                       detail::is_type_specialisation_of_v<remove_cvref_t<decltype(std::get<1>(mate))>, std::optional>) &&
                      std::Integral<remove_cvref_t<decltype(std::get<2>(mate))>>),
                      "The mate object must be a std::tuple of size 3 with "
                      "1) a std::ranges::ForwardRange with a value_type modelling seqan3::Alphabet, "
                      "2) a std::Integral or std::optional<std::Integral>, and "
                      "3) a std::Integral.");

        if constexpr (std::Integral<remove_cvref_t<decltype(std::get<0>(mate))>> ||
                      detail::is_type_specialisation_of_v<remove_cvref_t<decltype(std::get<0>(mate))>, std::optional>)
            static_assert(!detail::decays_to_ignore_v<header_type>,
                          "If you give indices as mate reference id information the header must also be present.");

        static_assert(std::Same<remove_cvref_t<tag_dict_type>, sam_tag_dictionary>,
                      "The tag_dict object must be of type seqan3::sam_tag_dictionary.");

        // ---------------------------------------------------------------------
        // logical Requirements
        // ---------------------------------------------------------------------
        if constexpr (!detail::decays_to_ignore_v<header_type> &&
                      !detail::decays_to_ignore_v<ref_id_type> &&
                      !std::Integral<std::remove_reference_t<ref_id_type>> &&
                      !detail::is_type_specialisation_of_v<std::remove_reference_t<ref_id_type>, std::optional>)
        {

            if (options.sam_require_header && !std::ranges::empty(ref_id))
            {
                auto id_it = header.ref_dict.end();

                if constexpr (std::ranges::ContiguousRange<decltype(ref_id)> &&
                              std::ranges::SizedRange<decltype(ref_id)> &&
                              ForwardingRange<decltype(ref_id)>)
                {
                    id_it = header.ref_dict.find(std::span{std::ranges::data(ref_id), std::ranges::size(ref_id)});
                }
                else
                {
                    using header_ref_id_type = std::remove_reference_t<decltype(header.ref_ids()[0])>;

                    static_assert(ImplicitlyConvertibleTo<ref_id_type, header_ref_id_type>,
                                  "The ref_id type is not convertible to the reference id information stored in the "
                                  "reference dictionary of the header object.");

                    id_it = header.ref_dict.find(ref_id);
                }

                if (id_it == header.ref_dict.end()) // no reference id matched
                    throw format_error{detail::to_string("The ref_id '", ref_id, "' was not in the list of references:",
                                                         header.ref_ids())};
            }
        }

        if (ref_offset.has_value() && (ref_offset.value() + 1) < 0)
            throw format_error{"The ref_offset object must be an std::Integral >= 0."};

        // ---------------------------------------------------------------------
        // Writing the Header on first call
        // ---------------------------------------------------------------------
        if constexpr (!detail::decays_to_ignore_v<header_type>)
        {
            if (options.sam_require_header && !written_header)
            {
                write_header(stream, options, header);
                written_header = true;
            }
        }

        // ---------------------------------------------------------------------
        // Writing the Record
        // ---------------------------------------------------------------------
        seqan3::ostreambuf_iterator stream_it{stream};
        char const separator{'\t'};

        write_range(stream_it, std::forward<id_type>(id));

        stream << separator;

        stream << flag << separator;

        if constexpr (!detail::decays_to_ignore_v<ref_id_type>)
        {
            if constexpr (std::Integral<std::remove_reference_t<ref_id_type>>)
            {
                write_range(stream_it, (header.ref_ids())[ref_id]);
            }
            else if constexpr (detail::is_type_specialisation_of_v<std::remove_reference_t<ref_id_type>, std::optional>)
            {
                if (ref_id.has_value())
                    write_range(stream_it, (header.ref_ids())[ref_id.value()]);
                else
                    stream << '*';
            }
            else
            {
                write_range(stream_it, std::forward<ref_id_type>(ref_id));
            }
        }
        else
        {
            stream << '*';
        }

        stream << separator;

        // SAM is 1 based, 0 indicates unmapped read if optional is not set
        stream << (ref_offset.value_or(-1) + 1) << separator;

        stream << static_cast<unsigned>(mapq) << separator;

        if (!std::ranges::empty(get<0>(align)) && !std::ranges::empty(get<1>(align)))
        {
            // compute possible distance from alignment end to sequence end
            // which indicates soft clipping at the end.
            // This should be replace by a free count_gaps function for
            // aligned sequences which is more efficient if possible.
            size_t off_end{std::ranges::size(seq) - offset};
            for (auto chr : get<1>(align))
                if (chr == gap{})
                    ++off_end;
            off_end -= std::ranges::size(get<1>(align));

            write_range(stream_it, detail::get_cigar_string(std::forward<align_type>(align), offset, off_end));
        }
        else
        {
            stream << '*';
        }

        stream << separator;

        if constexpr (std::Integral<std::remove_reference_t<decltype(get<0>(mate))>>)
        {
            write_range(stream_it, (header.ref_ids())[get<0>(mate)]);
        }
        else if constexpr (detail::is_type_specialisation_of_v<std::remove_reference_t<decltype(get<0>(mate))>, std::optional>)
        {
            if (get<0>(mate).has_value())
                // value_or(0) instead of value() (which is equivalent here) as a
                // workaround for a ubsan false-positive in GCC8: https://gcc.gnu.org/bugzilla/show_bug.cgi?id=90058
                write_range(stream_it, header.ref_ids()[get<0>(mate).value_or(0)]);
            else
                stream << '*';
        }
        else
        {
            write_range(stream_it, get<0>(mate));
        }

        stream << separator;

        if constexpr (detail::is_type_specialisation_of_v<remove_cvref_t<decltype(get<1>(mate))>, std::optional>)
        {
            // SAM is 1 based, 0 indicates unmapped read if optional is not set
            stream << (get<1>(mate).value_or(-1) + 1) << separator;
        }
        else
        {
            stream << get<1>(mate) << separator;
        }

        stream << get<2>(mate) << separator;

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
    //!\brief A variable that tracks whether the content of header has been written or not.
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
        if (std::ranges::empty(field_value))
            stream_it = '*';
        else
            std::ranges::copy(field_value | view::to_char, stream_it);
    }

    /*!\brief Writes a field value to the stream.
     * \tparam stream_it_t The stream iterator type.
     *
     * \param[in,out] stream_it   The stream iterator to print to.
     * \param[in]     field_value The value to print.
     */
    template <typename stream_it_t>
    void write_range(stream_it_t & stream_it, char const * const field_value)
    {
        write_range(stream_it, std::string_view{field_value});
    }

    /*!\brief Writes a field value to the stream.
     * \tparam stream_t           The stream type.
     * \param[in,out] stream      The stream to print to.
     * \param[in]     field_value The value to print.
     */
    template <typename stream_t, Arithmetic field_type>
    void write_field(stream_t & stream, field_type field_value)
    {
        // TODO: replace this with to_chars for efficiency
        if constexpr (std::Same<field_type, int8_t> || std::Same<field_type, uint8_t>)
            stream << static_cast<int16_t>(field_value);
        else
            stream << field_value;
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
        auto stream_variant_fn = [this, &stream] (auto && arg) // helper to print an std::variant
        {
            using T = remove_cvref_t<decltype(arg)>;

            if constexpr (!Container<T> || std::Same<T, std::string>)
            {
                stream << arg;
            }
            else
            {
                if (arg.begin() != arg.end())
                {
                    for (auto it = arg.begin(); it != (arg.end() - 1); ++it)
                    {
                        write_field(stream, *it);
                        stream << ',';
                    }

                    write_field(stream, *(arg.end() - 1)); // write last value without trailing ','
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

    /*!\brief Writes the SAM header.
     * \tparam stream_t   The stream type.
     *
     * \param[in,out] stream  The stream to print to.
     * \param[in]     options The options to alter printing.
     * \param[in]     header  The header to print.
     *
     * \throws seqan3::format_error if the header object contains the wrong
     *         information or the contents are ill-formed.
     *
     * \details
     *
     * Before writing the header, the contents are checked for correctness
     * according to the rules of the official
     * [SAM format specifications](https://samtools.github.io/hts-specs/SAMv1.pdf).
     */
    template <typename stream_t, typename ref_ids_type>
    void write_header(stream_t & stream,
                      alignment_file_output_options const & options,
                      alignment_file_header<ref_ids_type> & header)
    {
        // -----------------------------------------------------------------
        // Check Header
        // -----------------------------------------------------------------

        // (@HD) Check header line
        // The format version string will be taken from the local member variable
        if (!header.sorting.empty() &&
            !(header.sorting == "unknown"   ||
              header.sorting == "unsorted"  ||
              header.sorting == "queryname" ||
              header.sorting == "coordinate" ))
            throw format_error{"SAM format error: The header.sorting member must be "
                               "one of [unknown, unsorted, queryname, coordinate]."};

        if (!header.grouping.empty() &&
            !(header.grouping == "none"   ||
              header.grouping == "query"  ||
              header.grouping == "reference"))
            throw format_error{"SAM format error: The header.grouping member must be "
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
        seqan3::ostreambuf_iterator stream_it{stream};

        // (@HD) Write header line [required].
        stream << "@HD\tVN:";
        stream << format_sam::format_version;

        if (!header.sorting.empty())
            stream << "\tSO:" << header.sorting;

        if (!header.subsorting.empty())
            stream << "\tSS:" << header.subsorting;

        if (!header.grouping.empty())
            stream << "\tGO:" << header.grouping;

        detail::write_eol(stream_it, options.add_carriage_return);

        // (@SQ) Write Reference Sequence Dictionary lines [required].
        for (auto const & [ref_name, ref_info] : std::view::zip(header.ref_ids(), header.ref_id_info))
        {
            stream << "@SQ\tSN:";

            std::ranges::copy(ref_name, stream_it);

            stream << "\tLN:" << get<0>(ref_info);

            if (!get<1>(ref_info).empty())
                stream << "\t" << get<1>(ref_info);

            detail::write_eol(stream_it, options.add_carriage_return);
        }

        // Write read group (@RG) lines if specified.
        for (auto const & read_group : header.read_groups)
        {
            stream << "@RG"
                   << "\tID:" << get<0>(read_group);

            if (!get<1>(read_group).empty())
                stream << "\t" << get<1>(read_group);

            detail::write_eol(stream_it, options.add_carriage_return);
        }

        // Write program (@PG) lines if specified.
        for (auto const & program : header.program_infos)
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
        for (auto const & comment : header.comments)
        {
            stream << "@CO\t" << comment;
            detail::write_eol(stream_it, options.add_carriage_return);
        }
    }
};

} // namespace seqan3::detail
