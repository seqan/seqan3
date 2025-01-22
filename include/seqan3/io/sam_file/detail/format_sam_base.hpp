// SPDX-FileCopyrightText: 2006-2025 Knut Reinert & Freie Universität Berlin
// SPDX-FileCopyrightText: 2016-2025 Knut Reinert & MPI für molekulare Genetik
// SPDX-License-Identifier: BSD-3-Clause

/*!\file
 * \brief Provides the seqan3::format_sam_base that can be inherited from.
 * \author Svenja Mehringer <svenja.mehringer AT fu-berlin.de>
 * \author Lydia Buntrock <lydia.buntrock AT fu-berlin.de>
 */

#pragma once

#include <climits>
#include <ranges>
#include <string>
#include <vector>

#include <seqan3/alphabet/views/char_to.hpp>
#include <seqan3/core/debug_stream/range.hpp>
#include <seqan3/core/range/detail/misc.hpp>
#include <seqan3/io/detail/misc.hpp>
#include <seqan3/io/sam_file/detail/cigar.hpp>
#include <seqan3/io/sam_file/header.hpp>
#include <seqan3/io/sam_file/output_format_concept.hpp>
#include <seqan3/utility/detail/type_name_as_string.hpp>
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
    format_sam_base() = default;                                    //!< Defaulted.
    format_sam_base(format_sam_base const &) = default;             //!< Defaulted.
    format_sam_base & operator=(format_sam_base const &) = default; //!< Defaulted.
    format_sam_base(format_sam_base &&) = default;                  //!< Defaulted.
    format_sam_base & operator=(format_sam_base &&) = default;      //!< Defaulted.
    ~format_sam_base() = default;                                   //!< Defaulted.

    //!\}

    //!\brief The format version string.
    static constexpr std::array format_version{'1', '.', '6'};

    //!\brief A buffer used when parsing arithmetic values with std::from_chars.
    std::array<char, 316> arithmetic_buffer{}; // Doubles can be up to 316 characters

    //!\brief A variable that tracks whether the content of header has been written or not.
    bool header_was_written{false};

    //!\brief Tracks whether reference information (\@SQ tag) were found in the SAM header
    bool ref_info_present_in_header{false};

    template <typename ref_id_type, typename ref_id_tmp_type, typename header_type, typename ref_seqs_type>
    void check_and_assign_ref_id(ref_id_type & ref_id,
                                 ref_id_tmp_type & ref_id_tmp,
                                 header_type & header,
                                 ref_seqs_type & /*tag*/);

    int32_t soft_clipping_at_front(std::vector<cigar> const & cigar_vector) const;

    template <typename stream_view_type, std::ranges::forward_range target_range_type>
    void read_forward_range_field(stream_view_type && stream_view, target_range_type & target);

    template <std::ranges::forward_range target_range_type>
    void read_forward_range_field(std::string_view const str, target_range_type & target);

    template <arithmetic arithmetic_target_type>
    void read_arithmetic_field(std::string_view const & str, arithmetic_target_type & arithmetic_target);

    template <typename stream_view_type, typename ref_ids_type, typename ref_seqs_type>
    void read_header(stream_view_type && stream_view,
                     sam_file_header<ref_ids_type> & hdr,
                     ref_seqs_type & /*ref_id_to_pos_map*/);

    template <typename stream_t, typename header_type>
    void write_header(stream_t & stream, sam_file_output_options const & options, header_type & header);
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
template <typename ref_id_type, typename ref_id_tmp_type, typename header_type, typename ref_seqs_type>
inline void format_sam_base::check_and_assign_ref_id(ref_id_type & ref_id,
                                                     ref_id_tmp_type & ref_id_tmp,
                                                     header_type & header,
                                                     ref_seqs_type & /*tag*/)
{
    if (!std::ranges::empty(ref_id_tmp)) // otherwise the std::optional will not be filled
    {
        auto search = header.ref_dict.find(ref_id_tmp);

        if (search == header.ref_dict.end())
        {
            if constexpr (detail::decays_to_ignore_v<ref_seqs_type>) // no reference information given
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

/*!\brief Returns the soft clipping value at the front of the \p cigar_vector or 0 if none present.
 * \param[in] cigar_vector The cigar information to parse for soft-clipping.
 */
inline int32_t format_sam_base::soft_clipping_at_front(std::vector<cigar> const & cigar_vector) const
{
    int32_t sc_front{};

    // Checks if the given index in the cigar vector is a soft clip.
    auto soft_clipping_at = [&](size_t const index)
    {
        return cigar_vector[index] == 'S'_cigar_operation;
    };
    // Checks if the given index in the cigar vector is a hard clip.
    auto hard_clipping_at = [&](size_t const index)
    {
        return cigar_vector[index] == 'H'_cigar_operation;
    };
    // Checks if the given cigar vector as at least min_size many elements.
    auto vector_size_at_least = [&](size_t const min_size)
    {
        return cigar_vector.size() >= min_size;
    };
    // Returns the cigar count of the ith cigar element in the given cigar vector.
    auto cigar_count_at = [&](size_t const index)
    {
        return get<0>(cigar_vector[index]);
    };

    // check for soft clipping at the first two positions
    if (vector_size_at_least(1) && soft_clipping_at(0))
        sc_front = cigar_count_at(0);
    else if (vector_size_at_least(2) && hard_clipping_at(0) && soft_clipping_at(1))
        sc_front = cigar_count_at(1);

    return sc_front;
}

/*!\brief Reads a range by copying from stream_view to target, converting values with seqan3::views::char_to.
 * \tparam stream_view_type  The type of the stream as a view.
 * \tparam target_range_type The type of range to parse from input; must model std::ranges::forward_range.
 *
 * \param[in, out] stream_view  The stream view to iterate over.
 * \param[out]     target       The range to store the parsed sequence.
 */
template <typename stream_view_type, std::ranges::forward_range target_range_type>
inline void format_sam_base::read_forward_range_field(stream_view_type && stream_view, target_range_type & target)
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
                              std::back_inserter(target));
        }
    }
}

/*!\brief Reads from `str` to `target`, converting values with seqan3::views::char_to.
 * \tparam target_range_type The type of range to parse from input; must model std::ranges::forward_range.
 *
 * \param[in, out] str          The string_view to parse.
 * \param[out]     target       The range to store the parsed sequence.
 */
template <std::ranges::forward_range target_range_type>
inline void format_sam_base::read_forward_range_field(std::string_view const str, target_range_type & target)
{
    if (str.size() == 1 && str[0] == '*') // '*' denotes empty field
        return;

    if constexpr (std::assignable_from<target_range_type, std::string_view>)
    {
        target = str;
    }
    else
    {
        target.resize(str.size());
        for (size_t i = 0; i < str.size(); ++i)
            target[i] = assign_char_to(str[i], std::ranges::range_value_t<target_range_type>{});
    }
}

/*!\brief Reads arithmetic fields using std::from_chars.
 * \tparam arithmetic_target_type      The type of value to parse from input; must model seqan3::arithmetic.
 *
 * \param[in, out] str  The string_view to parse.
 * \param[out] arithmetic_target The arithmetic value object to store the parsed value.
 *
 * \throws seqan3::format_error if the character sequence in str cannot be successfully converted to a value
 *         of type arithmetic_target_type.
 */
template <arithmetic arithmetic_target_type>
inline void format_sam_base::read_arithmetic_field(std::string_view const & str,
                                                   arithmetic_target_type & arithmetic_target)
{
    std::from_chars_result res = std::from_chars(str.begin(), str.end(), arithmetic_target);

    if (res.ec == std::errc::invalid_argument || res.ptr != str.end())
        throw format_error{std::string("[CORRUPTED SAM FILE] The string '") + std::string(str.begin(), str.end())
                           + "' could not be cast into type " + detail::type_name_as_string<arithmetic_target_type>};

    if (res.ec == std::errc::result_out_of_range)
        throw format_error{std::string("[CORRUPTED SAM FILE] Casting '") + std::string(str.begin(), str.end())
                           + "' into type " + detail::type_name_as_string<arithmetic_target_type>
                           + " would cause an overflow."};
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
 * The function throws a seqan3::format_error if the format is not in a correct state (e.g. required fields are not
 * given), but throwing might occur downstream of the actual error.
 *
 * Any user-defined tags are not checked for correctness ([TAG]:[VALUE]) and are stored as strings:
 * * HD: seqan3::sam_file_header::user_tags
 * * SQ: seqan3::sam_file_header::ref_id_info
 * * RG: seqan3::sam_file_header::read_groups
 * * PG: seqan3::sam_file_header::program_infos / seqan3::sam_file_program_info_t::user_tags
 */
template <typename stream_view_type, typename ref_ids_type, typename ref_seqs_type>
inline void format_sam_base::read_header(stream_view_type && stream_view,
                                         sam_file_header<ref_ids_type> & hdr,
                                         ref_seqs_type & /*ref_id_to_pos_map*/)
{
    auto it = std::ranges::begin(stream_view);
    auto end = std::ranges::end(stream_view);
    std::string string_buffer{};

    auto make_tag = [](uint8_t char1, uint8_t char2) constexpr
    {
        return static_cast<uint16_t>(char1) | (static_cast<uint16_t>(char2) << CHAR_BIT);
    };

    std::array<char, 2> raw_tag{};

    auto parse_and_make_tag = [&]()
    {
        raw_tag[0] = *it;
        ++it;
        raw_tag[1] = *it;
        ++it;
        return make_tag(raw_tag[0], raw_tag[1]);
    };

    auto take_until_predicate = [&it, &string_buffer](auto const & predicate)
    {
        string_buffer.clear();
        while (!predicate(*it))
        {
            string_buffer.push_back(*it);
            ++it;
        }
    };

    auto skip_until_predicate = [&it](auto const & predicate)
    {
        while (!predicate(*it))
            ++it;
    };

    auto copy_next_tag_value_into_buffer = [&]()
    {
        skip_until_predicate(is_char<':'>);
        ++it; // skip :
        take_until_predicate(is_char<'\t'> || is_char<'\n'>);
    };

    // Some tags are not parsed individually. Instead, these are simply copied into a std::string.
    // Multiple tags must be separated by a `\t`, hence we prepend a tab to the string, except the first time.
    // Alternatively, we could always append a `\t`, but this would have the side effect that we might need to trim a
    // trailing tab after parsing all tags via `pop_back()`.
    // Unfortunately, we do not know when we are parsing the last tag (and in this case just not append a tab),
    // because even though we can check if the line ends in a `\n`, it is not guaranteed that the last tag of the
    // line is passed to this lambda. For example, the line might end with a tag that is properly parsed, such as `ID`.
    auto parse_and_append_unhandled_tag_to_string = [&](std::string & value, std::array<char, 2> raw_tag)
    {
        take_until_predicate(is_char<'\t'> || is_char<'\n'>);
        if (!value.empty())
            value.push_back('\t');
        value.push_back(raw_tag[0]);
        value.push_back(raw_tag[1]);
        read_forward_range_field(string_buffer, value);
    };

    while (it != end && is_char<'@'>(*it))
    {
        ++it; // skip @

        switch (parse_and_make_tag())
        {
        case make_tag('H', 'D'): // HD (header) tag
        {
            // All tags can appear in any order, VN is the only required tag
            while (is_char<'\t'>(*it))
            {
                ++it; // skip tab
                std::string * header_entry{nullptr};

                switch (parse_and_make_tag())
                {
                case make_tag('V', 'N'): // parse required VN (version) tag
                {
                    header_entry = std::addressof(hdr.format_version);
                    break;
                }
                case make_tag('S', 'O'): // SO (sorting) tag
                {
                    header_entry = std::addressof(hdr.sorting);
                    break;
                }
                case make_tag('S', 'S'): // SS (sub-order) tag
                {
                    header_entry = std::addressof(hdr.subsorting);
                    break;
                }
                case make_tag('G', 'O'): // GO (grouping) tag
                {
                    header_entry = std::addressof(hdr.grouping);
                    break;
                }
                default: // unknown/user tag
                {
                    parse_and_append_unhandled_tag_to_string(hdr.user_tags, raw_tag);
                }
                }

                if (header_entry != nullptr)
                {
                    copy_next_tag_value_into_buffer();
                    read_forward_range_field(string_buffer, *header_entry);
                }
            }
            ++it; // skip newline

            if (hdr.format_version.empty())
                throw format_error{std::string{"The required VN tag in @HD is missing."}};

            break;
        }

        case make_tag('S', 'Q'): // SQ (sequence dictionary) tag
        {
            ref_info_present_in_header = true;
            std::ranges::range_value_t<decltype(hdr.ref_ids())> id;
            std::optional<int32_t> sequence_length{};
            std::tuple<int32_t, std::string> info{};

            // All tags can appear in any order, SN and LN are required tags
            while (is_char<'\t'>(*it))
            {
                ++it; // skip tab

                switch (parse_and_make_tag())
                {
                case make_tag('S', 'N'): // parse required SN (sequence name) tag
                {
                    copy_next_tag_value_into_buffer();
                    read_forward_range_field(string_buffer, id);
                    break;
                }
                case make_tag('L', 'N'): // parse required LN (length) tag
                {
                    int32_t sequence_length_tmp{};
                    copy_next_tag_value_into_buffer();
                    read_arithmetic_field(string_buffer, sequence_length_tmp);
                    sequence_length = sequence_length_tmp;
                    break;
                }
                default: // Any other tag
                {
                    parse_and_append_unhandled_tag_to_string(get<1>(info), raw_tag);
                }
                }
            }
            ++it; // skip newline

            if (id.empty())
                throw format_error{std::string{"The required SN tag in @SQ is missing."}};
            if (!sequence_length.has_value())
                throw format_error{std::string{"The required LN tag in @SQ is missing."}};
            if (sequence_length.value() <= 0)
                throw format_error{std::string{"The value of LN in @SQ must be positive."}};

            get<0>(info) = sequence_length.value();
            // If reference information was given, the ids exist and we can fill ref_dict directly.
            // If not, we need to update the ids first and fill the reference dictionary afterwards.
            if constexpr (!detail::decays_to_ignore_v<ref_seqs_type>) // reference information given
            {
                auto id_it = hdr.ref_dict.find(id);

                if (id_it == hdr.ref_dict.end())
                    throw format_error{detail::to_string("Unknown reference name '",
                                                         id,
                                                         "' found in SAM header ",
                                                         "(header.ref_ids(): ",
                                                         hdr.ref_ids(),
                                                         ").")};

                auto & given_ref_info = hdr.ref_id_info[id_it->second];

                if (std::get<0>(given_ref_info) != std::get<0>(info))
                    throw format_error{"Provided and header-based reference length differ."};

                hdr.ref_id_info[id_it->second] = std::move(info);
            }
            else
            {
                static_assert(!detail::is_type_specialisation_of_v<decltype(hdr.ref_ids()), std::deque>,
                              "The range over reference ids must be of type std::deque such that pointers are not "
                              "invalidated.");

                hdr.ref_ids().push_back(id);
                hdr.ref_id_info.push_back(info);
                hdr.ref_dict[(hdr.ref_ids())[(hdr.ref_ids()).size() - 1]] = (hdr.ref_ids()).size() - 1;
            }
            break;
        }

        case make_tag('R', 'G'): // RG (read group) tag
        {
            std::pair<std::string, std::string> tmp{};

            // All tags can appear in any order, SN and LN are required tags
            while (is_char<'\t'>(*it))
            {
                ++it; // skip tab

                switch (parse_and_make_tag())
                {
                case make_tag('I', 'D'): // parse required ID tag
                {
                    copy_next_tag_value_into_buffer();
                    read_forward_range_field(string_buffer, get<0>(tmp));
                    break;
                }
                default: // Any other tag
                {
                    parse_and_append_unhandled_tag_to_string(get<1>(tmp), raw_tag);
                }
                }
            }
            ++it; // skip newline

            if (get<0>(tmp).empty())
                throw format_error{std::string{"The required ID tag in @RG is missing."}};

            hdr.read_groups.emplace_back(std::move(tmp));
            break;
        }

        case make_tag('P', 'G'): // PG (program) tag
        {
            typename sam_file_header<ref_ids_type>::program_info_t tmp{};

            // All tags can appear in any order, ID is the only required tag
            while (is_char<'\t'>(*it))
            {
                ++it; // skip tab
                std::string * program_info_entry{nullptr};

                switch (parse_and_make_tag())
                {
                case make_tag('I', 'D'): // read required ID tag
                {
                    program_info_entry = std::addressof(tmp.id);
                    break;
                }
                case make_tag('P', 'N'): // PN (program name) tag
                {
                    program_info_entry = std::addressof(tmp.name);
                    break;
                }
                case make_tag('P', 'P'): // PP (previous program) tag
                {
                    program_info_entry = std::addressof(tmp.previous);
                    break;
                }
                case make_tag('C', 'L'): // CL (command line) tag
                {
                    program_info_entry = std::addressof(tmp.command_line_call);
                    break;
                }
                case make_tag('D', 'S'): // DS (description) tag
                {
                    program_info_entry = std::addressof(tmp.description);
                    break;
                }
                case make_tag('V', 'N'): // VN (version) tag
                {
                    program_info_entry = std::addressof(tmp.version);
                    break;
                }
                default: // unsupported header tag
                {
                    parse_and_append_unhandled_tag_to_string(tmp.user_tags, raw_tag);
                }
                }

                if (program_info_entry != nullptr)
                {
                    copy_next_tag_value_into_buffer();
                    read_forward_range_field(string_buffer, *program_info_entry);
                }
            }
            ++it; // skip newline

            if (tmp.id.empty())
                throw format_error{std::string{"The required ID tag in @PG is missing."}};

            hdr.program_infos.emplace_back(std::move(tmp));
            break;
        }

        case make_tag('C', 'O'): // CO (comment) tag
        {
            ++it; // skip tab
            std::string tmp;
            take_until_predicate(is_char<'\n'>);
            read_forward_range_field(string_buffer, tmp);
            ++it; // skip newline
            hdr.comments.emplace_back(std::move(tmp));
            break;
        }

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
template <typename stream_t, typename header_type>
inline void
format_sam_base::write_header(stream_t & stream, sam_file_output_options const & options, header_type & header)
{
    if constexpr (!detail::decays_to_ignore_v<header_type>)
    {
        // -----------------------------------------------------------------
        // Check Header
        // -----------------------------------------------------------------

        // (@HD) Check header line
        // The format version string will be taken from the local member variable
        if (!header.sorting.empty()
            && !(header.sorting == "unknown" || header.sorting == "unsorted" || header.sorting == "queryname"
                 || header.sorting == "coordinate"))
            throw format_error{"SAM format error: The header.sorting member must be "
                               "one of [unknown, unsorted, queryname, coordinate]."};

        if (!header.grouping.empty()
            && !(header.grouping == "none" || header.grouping == "query" || header.grouping == "reference"))
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
        std::ostreambuf_iterator stream_it{stream};

        // (@HD) Write header line [required].
        stream << "@HD\tVN:";
        std::ranges::copy(format_version, stream_it);

        if (!header.sorting.empty())
            stream << "\tSO:" << header.sorting;

        if (!header.subsorting.empty())
            stream << "\tSS:" << header.subsorting;

        if (!header.grouping.empty())
            stream << "\tGO:" << header.grouping;

        if (!header.user_tags.empty())
            stream << '\t' << header.user_tags;

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

            if (!program.user_tags.empty())
                stream << '\t' << program.user_tags;

            detail::write_eol(stream_it, options.add_carriage_return);
        }

        // Write comment (@CO) lines if specified.
        for (auto const & comment : header.comments)
        {
            stream << "@CO\t" << comment;
            detail::write_eol(stream_it, options.add_carriage_return);
        }
    }
}

} // namespace seqan3::detail
