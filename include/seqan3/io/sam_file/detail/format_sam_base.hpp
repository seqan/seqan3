// -----------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2020, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2020, Knut Reinert & MPI für molekulare Genetik
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at: https://github.com/seqan/seqan3/blob/master/LICENSE.md
// -----------------------------------------------------------------------------------------------------

/*!\file
 * \brief Provides the seqan3::format_sam_base that can be inherited from.
 * \author Svenja Mehringer <svenja.mehringer AT fu-berlin.de>
 */

#pragma once

#include <seqan3/std/algorithm>
#include <seqan3/std/charconv>
#include <seqan3/std/concepts>
#include <iterator>
#include <seqan3/std/ranges>
#include <string>
#include <vector>

#include <seqan3/alphabet/views/char_to.hpp>
#include <seqan3/alphabet/views/to_char.hpp>
#include <seqan3/core/debug_stream/detail/to_string.hpp>
#include <seqan3/core/detail/debug_stream_range.hpp>
#include <seqan3/core/detail/template_inspection.hpp>
#include <seqan3/core/range/detail/misc.hpp>
#include <seqan3/core/range/type_traits.hpp>
#include <seqan3/io/detail/ignore_output_iterator.hpp>
#include <seqan3/io/detail/istreambuf_view.hpp>
#include <seqan3/io/detail/misc.hpp>
#include <seqan3/io/detail/take_until_view.hpp>
#include <seqan3/io/sam_file/detail/cigar.hpp>
#include <seqan3/io/sam_file/header.hpp>
#include <seqan3/io/sam_file/input_format_concept.hpp>
#include <seqan3/io/sam_file/input_options.hpp>
#include <seqan3/io/sam_file/output_format_concept.hpp>
#include <seqan3/io/sam_file/output_options.hpp>
#include <seqan3/io/sam_file/sam_tag_dictionary.hpp>
#include <seqan3/io/sequence_file/output_options.hpp>
#include <seqan3/utility/detail/exposition_only_concept.hpp>
#include <seqan3/utility/detail/type_name_as_string.hpp>
#include <seqan3/utility/tuple/concept.hpp>
#include <seqan3/utility/views/repeat_n.hpp>
#include <seqan3/utility/views/slice.hpp>
#include <seqan3/utility/views/zip.hpp>

namespace seqan3::detail
{

/*!\brief The alignment base format.
 * \ingroup io_sam_file
 *
 * \details
 *
 * Since the SAM and BAM format share a lot of functionality, this abstract base class defines common member variables
 * and functions that are used in both formats.
 */
class format_sam_base
{
protected:
    /*!\name Constructors, destructor and assignment
     * \{
     */
    format_sam_base() = default; //!< Defaulted.
    format_sam_base(format_sam_base const &) = default; //!< Defaulted.
    format_sam_base & operator=(format_sam_base const &) = default; //!< Defaulted.
    format_sam_base(format_sam_base &&) = default; //!< Defaulted.
    format_sam_base & operator=(format_sam_base &&) = default; //!< Defaulted.
    ~format_sam_base() = default; //!< Defaulted.

    //!\}

    //!\brief The format version string.
    static constexpr std::array format_version{'1', '.', '6'};

    //!\brief A buffer used when parsing arithmetic values with std::from_chars.
    std::array<char, 316> arithmetic_buffer{}; // Doubles can be up to 316 characters

    //!\brief A variable that tracks whether the content of header has been written or not.
    bool header_was_written{false};

    //!\brief Tracks whether reference information (\@SR tag) were found in the SAM header
    bool ref_info_present_in_header{false};

    template <typename ref_id_type,
              typename ref_id_tmp_type,
              typename header_type,
              typename ref_seqs_type>
    void check_and_assign_ref_id(ref_id_type      & ref_id,
                                 ref_id_tmp_type  & ref_id_tmp,
                                 header_type      & header,
                                 ref_seqs_type    & /*tag*/);

    template <typename align_type, typename ref_seqs_type>
    void construct_alignment(align_type                           & align,
                             std::vector<cigar>                   & cigar_vector,
                             [[maybe_unused]] int32_t               rid,
                             [[maybe_unused]] ref_seqs_type       & ref_seqs,
                             [[maybe_unused]] int32_t               ref_start,
                             size_t                                 ref_length);

    void transfer_soft_clipping_to(std::vector<cigar> const & cigar_vector, int32_t & sc_begin, int32_t & sc_end) const;

    template <typename stream_view_type>
    void read_field(stream_view_type && stream_view, detail::ignore_t const & SEQAN3_DOXYGEN_ONLY(target));

    template <typename stream_view_t>
    void read_field(stream_view_t && stream_view, std::byte & byte_target);

    template <typename stream_view_type, std::ranges::forward_range target_range_type>
    void read_field(stream_view_type && stream_view, target_range_type & target);

    template <typename stream_view_t, arithmetic arithmetic_target_type>
    void read_field(stream_view_t && stream_view, arithmetic_target_type & arithmetic_target);

    template <typename stream_view_type, typename optional_value_type>
    void read_field(stream_view_type && stream_view, std::optional<optional_value_type> & target);

    template <typename stream_view_type, typename ref_ids_type, typename ref_seqs_type>
    void read_header(stream_view_type && stream_view,
                     sam_file_header<ref_ids_type> & hdr,
                     ref_seqs_type & /*ref_id_to_pos_map*/);

    template <typename stream_t, typename ref_ids_type>
    void write_header(stream_t & stream,
                      sam_file_output_options const & options,
                      sam_file_header<ref_ids_type> & header);
};

/*!\brief Checks for known reference ids or adds a new reference is and assigns a reference id to `ref_id`.
 * \tparam ref_id_type         The type of the reference id (usually a views::type_reduce over ref_id_tmp_type).
 * \tparam ref_id_tmp_type     The type of the temporary parsed id (same_as type as reference ids in header).
 * \tparam header_type         The type of the alignment header.
 * \tparam ref_seqs_type       A tag whether the reference information were given or not (std::ignore or not).
 *
 * \param[out]     ref_id      The reference id to be filled.
 * \param[in]      ref_id_tmp  The temporary of the parsed reference id.
 * \param[in, out] header      The header object that stores the reference id information.
 */
template <typename ref_id_type,
          typename ref_id_tmp_type,
          typename header_type,
          typename ref_seqs_type>
inline void format_sam_base::check_and_assign_ref_id(ref_id_type      & ref_id,
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

/*!\brief Transfer soft clipping information from the \p cigar_vector to \p sc_begin and \p sc_end.
 * \param[in] cigar_vector The cigar information to parse for soft-clipping.
 * \param[out] sc_begin    The soft clipping at the beginning of the alignment to set.
 * \param[out] sc_end      The soft clipping at the end of the alignment to set.
 */
inline void format_sam_base::transfer_soft_clipping_to(std::vector<cigar> const & cigar_vector,
                                                       int32_t & sc_begin,
                                                       int32_t & sc_end) const
{
    // Checks if the given index in the cigar vector is a soft clip.
    auto soft_clipping_at = [&] (size_t const index) { return cigar_vector[index] == 'S'_cigar_operation; };
    // Checks if the given index in the cigar vector is a hard clip.
    auto hard_clipping_at = [&] (size_t const index) { return cigar_vector[index] == 'H'_cigar_operation; };
    // Checks if the given cigar vector as at least min_size many elements.
    auto vector_size_at_least = [&] (size_t const min_size) { return cigar_vector.size() >= min_size; };
    // Returns the cigar count of the ith cigar element in the given cigar vector.
    auto cigar_count_at = [&] (size_t const index) { return get<0>(cigar_vector[index]); };

    // check for soft clipping at the first two positions
    if (vector_size_at_least(1) && soft_clipping_at(0))
        sc_begin = cigar_count_at(0);
    else if (vector_size_at_least(2) && hard_clipping_at(0) && soft_clipping_at(1))
        sc_begin = cigar_count_at(1);

    // Check for soft clipping at the last two positions. But only if the vector size has at least 2, respectively
    // 3 elements. Accordingly, if the following arithmetics overflow they are protected by the corresponding
    // if expressions below.
    auto last_index = cigar_vector.size() - 1;
    auto second_last_index = last_index - 1;

    if (vector_size_at_least(2) && soft_clipping_at(last_index))
        sc_end = cigar_count_at(last_index);
    else if (vector_size_at_least(3) && hard_clipping_at(last_index) && soft_clipping_at(second_last_index))
        sc_end = cigar_count_at(second_last_index);
}

/*!\brief Construct the field::alignment depending on the given information.
 * \tparam align_type      The alignment type.
 * \tparam ref_seqs_type   The type of reference sequences (might decay to ignore).
 * \param[in,out] align    The alignment (pair of aligned sequences) to fill.
 * \param[in] cigar_vector The cigar information to convert to an alignment.
 * \param[in] rid          The index of the reference sequence in header.ref_ids().
 * \param[in] ref_seqs     The reference sequence information.
 * \param[in] ref_start    The start position of the alignment in the reference sequence.
 * \param[in] ref_length   The length of the aligned reference sequence.
 */
template <typename align_type, typename ref_seqs_type>
inline void format_sam_base::construct_alignment(align_type                           & align,
                                                 std::vector<cigar>                   & cigar_vector,
                                                 [[maybe_unused]] int32_t               rid,
                                                 [[maybe_unused]] ref_seqs_type       & ref_seqs,
                                                 [[maybe_unused]] int32_t               ref_start,
                                                 size_t                                 ref_length)
{
    if (rid > -1 && ref_start > -1 &&       // read is mapped
        !cigar_vector.empty() &&                   // alignment field was not empty
        !std::ranges::empty(get<1>(align))) // seq field was not empty
    {
        if constexpr (!detail::decays_to_ignore_v<ref_seqs_type>)
        {
            assert(static_cast<size_t>(ref_start + ref_length) <= std::ranges::size(ref_seqs[rid]));
            // copy over unaligned reference sequence part
            assign_unaligned(get<0>(align), ref_seqs[rid] | views::slice(ref_start, ref_start + ref_length));
        }
        else
        {
            using unaligned_t = std::remove_cvref_t<detail::unaligned_seq_t<decltype(get<0>(align))>>;
            auto dummy_seq    = views::repeat_n(std::ranges::range_value_t<unaligned_t>{}, ref_length)
                              | std::views::transform(detail::access_restrictor_fn{});
            static_assert(std::same_as<unaligned_t, decltype(dummy_seq)>,
                          "No reference information was given so the type of the first alignment tuple position"
                          "must have an unaligned sequence type of a dummy sequence ("
                          "views::repeat_n(dna5{}, size_t{}) | "
                          "std::views::transform(detail::access_restrictor_fn{}))");

            assign_unaligned(get<0>(align), dummy_seq); // assign dummy sequence
        }

        // insert gaps according to the cigar information
        detail::alignment_from_cigar(align, cigar_vector);
    }
    else // not enough information for an alignment, assign an empty view/dummy_sequence
    {
        if constexpr (!detail::decays_to_ignore_v<ref_seqs_type>) // reference info given
        {
            assert(std::ranges::size(ref_seqs) > 0); // we assume that the given ref info is not empty
            assign_unaligned(get<0>(align), ref_seqs[0] | views::slice(0, 0));
        }
        else
        {
            using unaligned_t = std::remove_cvref_t<detail::unaligned_seq_t<decltype(get<0>(align))>>;
            assign_unaligned(get<0>(align), views::repeat_n(std::ranges::range_value_t<unaligned_t>{}, 0)
                                            | std::views::transform(detail::access_restrictor_fn{}));
        }
    }
}

/*!\brief Decays to detail::consume for std::ignore.
 * \tparam stream_view_type  The type of the stream as a view.
 *
 * \param[in, out] stream_view  The stream view to consume.
 * \param[in]      target       A std::ignore placeholder.
 */
template <typename stream_view_type>
inline void format_sam_base::read_field(stream_view_type && stream_view,
                                        detail::ignore_t const & SEQAN3_DOXYGEN_ONLY(target))
{
    detail::consume(stream_view);
}

/*!\brief Reads std::byte fields using std::from_chars.
 * \tparam stream_view_t The type of the stream as a view.
 *
 * \param[in, out] stream_view  The stream view to iterate over.
 * \param[out] byte_target The std::byte object to store the parsed value.
 *
 * \throws seqan3::format_error if the character sequence in stream_view cannot be successfully converted to a value
 *         of type std::byte.
 */
template <typename stream_view_t>
inline void format_sam_base::read_field(stream_view_t && stream_view, std::byte & byte_target)
{
    // unfortunately std::from_chars only accepts char const * so we need a buffer.
    auto [ignore, end] = std::ranges::copy(stream_view, arithmetic_buffer.data());
    (void) ignore;

    uint8_t byte{};
    // std::from_chars cannot directly parse into a std::byte
    std::from_chars_result res = std::from_chars(arithmetic_buffer.begin(), end, byte, 16);

    if (res.ec == std::errc::invalid_argument || res.ptr != end)
        throw format_error{std::string("[CORRUPTED SAM FILE] The string '") +
                                       std::string(arithmetic_buffer.begin(), end) +
                                       "' could not be cast into type uint8_t."};

    if (res.ec == std::errc::result_out_of_range)
        throw format_error{std::string("[CORRUPTED SAM FILE] Casting '") + std::string(arithmetic_buffer.begin(), end) +
                                       "' into type uint8_t would cause an overflow."};
    byte_target = std::byte{byte};
}

/*!\brief Reads a range by copying from stream_view to target, converting values with seqan3::views::char_to.
 * \tparam stream_view_type  The type of the stream as a view.
 * \tparam target_range_type The type of range to parse from input; must model std::ranges::forward_range.
 *
 * \param[in, out] stream_view  The stream view to iterate over.
 * \param[out]     target       The range to store the parsed sequence.
 */
template <typename stream_view_type, std::ranges::forward_range target_range_type>
inline void format_sam_base::read_field(stream_view_type && stream_view, target_range_type & target)
{
    using target_range_value_t = std::ranges::range_value_t<target_range_type>;
    using begin_iterator_t = std::ranges::iterator_t<stream_view_type>;
    using end_iterator_t = std::ranges::sentinel_t<stream_view_type>;

    // Note that we need to cache the begin iterator since the stream_view is an input range that may be consuming
    // and in that case might read `past-the-end` on a second call of std::ranges::begin.
    if (auto it = std::ranges::begin(stream_view); it != std::ranges::end(stream_view))
    {
        // Write to target if field does not represent an empty string, denoted as single '*' character.
        if (char c = *it; !(++it == std::ranges::end(stream_view) && c == '*'))
        {
            target.push_back(seqan3::assign_char_to(c, target_range_value_t{}));
            std::ranges::copy(std::ranges::subrange<begin_iterator_t, end_iterator_t>{it, std::ranges::end(stream_view)}
                              | views::char_to<target_range_value_t>,
                              std::cpp20::back_inserter(target));
        }
    }
}

/*!\brief Reads arithmetic fields using std::from_chars.
 * \tparam stream_view_t The type of the stream as a view.
 * \tparam arithmetic_target_type      The type of value to parse from input; must model seqan3::arithmetic.
 *
 * \param[in, out] stream_view  The stream view to iterate over.
 * \param[out] arithmetic_target The arithmetic value object to store the parsed value.
 *
 * \throws seqan3::format_error if the character sequence in stream_view cannot be successfully converted to a value
 *         of type arithmetic_target_type.
 */
template <typename stream_view_t, arithmetic arithmetic_target_type>
inline void format_sam_base::read_field(stream_view_t && stream_view, arithmetic_target_type & arithmetic_target)
{
    // unfortunately std::from_chars only accepts char const * so we need a buffer.
    auto [ignore, end] = std::ranges::copy(stream_view, arithmetic_buffer.data());
    (void) ignore;
    std::from_chars_result res = std::from_chars(arithmetic_buffer.begin(), end, arithmetic_target);

    if (res.ec == std::errc::invalid_argument || res.ptr != end)
        throw format_error{std::string("[CORRUPTED SAM FILE] The string '") +
                                       std::string(arithmetic_buffer.begin(), end) +
                                       "' could not be cast into type " +
                                       detail::type_name_as_string<arithmetic_target_type>};

    if (res.ec == std::errc::result_out_of_range)
        throw format_error{std::string("[CORRUPTED SAM FILE] Casting '") + std::string(arithmetic_buffer.begin(), end) +
                                       "' into type " + detail::type_name_as_string<arithmetic_target_type> +
                                       " would cause an overflow."};
}

/*!\brief Delegate parsing of std::optional types to parsing of the inner value type.
 * \tparam stream_view_type     The type of the stream as a view.
 * \tparam optional_value_type  The inner type of a the std::optional type of \p target.
 *
 * \param[in, out] stream_view The stream view to iterate over.
 * \param[out] target The std::optional object to store the parsed value.
 *
 * \throws seqan3::format_error if the character sequence in stream_view cannot be successfully converted to a value
 *         of type target_type.
 */
template <typename stream_view_type, typename optional_value_type>
inline void format_sam_base::read_field(stream_view_type && stream_view, std::optional<optional_value_type> & target)
{
    optional_value_type tmp;
    read_field(std::forward<stream_view_type>(stream_view), tmp);
    target = tmp;
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
inline void format_sam_base::read_header(stream_view_type && stream_view,
                                         sam_file_header<ref_ids_type> & hdr,
                                         ref_seqs_type & /*ref_id_to_pos_map*/)
{
    auto it = std::ranges::begin(stream_view);
    auto end = std::ranges::end(stream_view);
    std::vector<char> string_buffer{};

    auto take_until_predicate = [&it, &string_buffer] (auto const & predicate)
    {
        string_buffer.clear();
        while (!predicate(*it))
        {
            string_buffer.push_back(*it);
            ++it;
        }
    };

    auto skip_tag_name_and_colon = [&it] ()
    {
        while (!is_char<':'>(*it))
            ++it;
        ++it;
    };

    auto parse_tag_value = [&] (auto & value) // helper function to parse the next tag value
    {
        skip_tag_name_and_colon();
        take_until_predicate(is_char<'\t'> || is_char<'\n'>);
        read_field(string_buffer, value);
    };

    auto parse_until_newline = [&] (auto & value)
    {
        take_until_predicate(is_char<'\n'>);
        read_field(string_buffer, value);
    };

    while (it != end && is_char<'@'>(*it))
    {
        ++it;                                                           // skip @

        switch (*it)
        {
            case 'H': /* HD (header) tag */
            {
                parse_tag_value(hdr.format_version);                    // parse required VN (version) tag

                // The SO, SS and GO tag are optional and can appear in any order
                while (is_char<'\t'>(*it))
                {
                    ++it;                                               // skip tab
                    std::string * who = std::addressof(hdr.grouping);

                    if (is_char<'S'>(*it))
                    {
                        ++it;                                           // skip S

                        if (is_char<'O'>(*it))                          // SO (sorting) tag
                            who = std::addressof(hdr.sorting);
                        else if (is_char<'S'>(*it))                     // SS (sub-order) tag
                            who = std::addressof(hdr.subsorting);
                        else
                            throw format_error{std::string{"Illegal SAM header tag: S"} + *it};
                    }
                    else if (!is_char<'G'>(*it))                        // GO (grouping) tag
                    {
                        throw format_error{std::string{"Illegal SAM header tag in @HG starting with: "} + *it};
                    }

                    parse_tag_value(*who);
                }
                ++it;                                                   // skip newline
            } break;

            case 'S': /* SQ (sequence dictionary) tag */
            {
                ref_info_present_in_header = true;
                std::ranges::range_value_t<decltype(hdr.ref_ids())> id;
                std::tuple<int32_t, std::string> info{};

                parse_tag_value(id);                                    // parse required SN (sequence name) tag
                ++it;                                                   // skip tab or newline
                parse_tag_value(get<0>(info));                          // parse required LN (length) tag

                if (is_char<'\t'>(*it))                                 // read rest of the tags
                {
                    ++it;                                               // skip tab
                    parse_until_newline(get<1>(info));
                }
                ++it;                                                   // skip newline

                /* If reference information was given, the ids exist and we can fill ref_dict directly.
                 * If not, wee need to update the ids first and fill the reference dictionary afterwards. */
                if constexpr (!detail::decays_to_ignore_v<ref_seqs_type>) // reference information given
                {
                    auto id_it = hdr.ref_dict.find(id);

                    if (id_it == hdr.ref_dict.end())
                        throw format_error{detail::to_string("Unknown reference name '", id, "' found in SAM header ",
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
            } break;

            case 'R': /* RG (read group) tag */
            {
                std::pair<std::string, std::string> tmp{};

                parse_tag_value(get<0>(tmp));                           // read required ID tag

                if (is_char<'\t'>(*it))                                 // read rest of the tags
                {
                    ++it;
                    parse_until_newline(get<1>(tmp));
                }
                ++it;                                                   // skip newline

                hdr.read_groups.emplace_back(std::move(tmp));
            } break;

            case 'P': /* PG (program) tag */
            {
                typename sam_file_header<ref_ids_type>::program_info_t tmp{};

                parse_tag_value(tmp.id);                                // read required ID tag

                // The PN, CL, PP, DS, VN are optional tags and can be given in any order.
                while (is_char<'\t'>(*it))
                {
                    ++it;                                               // skip tab
                    std::string * who = &tmp.version;

                    switch (*it)
                    {
                        case 'P':
                            ++it;                                       // skip P

                            if (is_char<'N'>(*it))                      // PN (program name) tag
                                who = &tmp.name;
                            else                                        // PP (previous program) tag
                                who = &tmp.previous;
                            break;

                        case 'C': // CL (command line) tag
                            who = &tmp.command_line_call;
                            break;

                        case 'D': // DS (description) tag
                            who = &tmp.description;
                            break;

                        case 'V': // VN (version) tag
                            // already set as default above
                            break;

                        default:
                            throw format_error{std::string{"Illegal SAM header tag starting with: "} + *it};
                    }

                    parse_tag_value(*who);
                }
                ++it;                                                   // skip newline

                hdr.program_infos.emplace_back(std::move(tmp));
            } break;

            case 'C': /* CO (comment) tag */
            {
                std::string tmp;
                ++it;                                                   // skip C
                ++it;                                                   // skip O
                ++it;                                                   // skip :
                parse_until_newline(tmp);
                ++it;                                                   // skip newline

                hdr.comments.emplace_back(std::move(tmp));
            } break;

            default:
                throw format_error{std::string{"Illegal SAM header tag starting with:"} + *it};
        }
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
inline void format_sam_base::write_header(stream_t & stream,
                                          sam_file_output_options const & options,
                                          sam_file_header<ref_ids_type> & header)
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
    std::cpp20::ostreambuf_iterator stream_it{stream};

    // (@HD) Write header line [required].
    stream << "@HD\tVN:";
    std::ranges::copy(format_version, stream_it);

    if (!header.sorting.empty())
        stream << "\tSO:" << header.sorting;

    if (!header.subsorting.empty())
        stream << "\tSS:" << header.subsorting;

    if (!header.grouping.empty())
        stream << "\tGO:" << header.grouping;

    detail::write_eol(stream_it, options.add_carriage_return);

    // (@SQ) Write Reference Sequence Dictionary lines [required].
    for (auto const & [ref_name, ref_info] : views::zip(header.ref_ids(), header.ref_id_info))
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

} // namespace seqan3::detail
