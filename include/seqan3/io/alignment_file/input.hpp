// -----------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2020, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2020, Knut Reinert & MPI für molekulare Genetik
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at: https://github.com/seqan/seqan3/blob/master/LICENSE.md
// -----------------------------------------------------------------------------------------------------

/*!\file
 * \brief [DEPRECATED] Provides seqan3::sam_file_input and corresponding traits classes.
 * \deprecated This header will be removed in 3.1.0; Please \#include <seqan3/io/sam_file/input.hpp> instead.
 * \author Svenja Mehringer <svenja.mehringer AT fu-berlin.de>
 */

#pragma once

#include <seqan3/io/sam_file/input.hpp>

#include <seqan3/io/alignment_file/input_options.hpp>

namespace seqan3
{
//!\deprecated Use seqan3::sam_file_input_default_traits instead.
template <typename ref_sequences_t = ref_info_not_given, typename ref_ids_t = std::deque<std::string>>
using alignment_file_input_default_traits SEQAN3_DEPRECATED_310
    = sam_file_input_default_traits<ref_sequences_t, ref_ids_t>;

//!\deprecated Use seqan3::sam_file_input instead.
template <
    sam_file_input_traits traits_type_ = sam_file_input_default_traits<>,
    detail::fields_specialisation selected_field_ids_ = fields<field::seq,
                                                               field::id,
                                                               field::offset,
                                                               field::ref_seq,
                                                               field::ref_id,
                                                               field::ref_offset,
                                                               field::alignment,
                                                               field::cigar,
                                                               field::mapq,
                                                               field::qual,
                                                               field::flag,
                                                               field::mate,
                                                               field::tags,
                                                               field::evalue,
                                                               field::bit_score,
                                                               field::header_ptr>,
    detail::type_list_of_sam_file_input_formats valid_formats_ = type_list<format_sam, format_bam>>
#if defined(__cpp_deduction_guides) && __cpp_deduction_guides >= 201907L
using alignment_file_input SEQAN3_DEPRECATED_310 = sam_file_input<traits_type_, selected_field_ids_, valid_formats_>;
#else // defined(__cpp_deduction_guides) && __cpp_deduction_guides >= 201907L
// This is a workaround for pre-c++20 as CTAD don't work on type aliases.
struct SEQAN3_DEPRECATED_310 alignment_file_input :
    public sam_file_input<traits_type_, selected_field_ids_, valid_formats_>
{
//!\cond
    using sam_file_input_t = sam_file_input<traits_type_, selected_field_ids_, valid_formats_>;
    using sam_file_input_t::sam_file_input_t;
    using sam_file_input_t::operator=;
//!\endcond
};

//!\cond
template <typename ...args_t>
alignment_file_input(args_t && ...args) ->
    alignment_file_input<typename decltype(sam_file_input{std::forward<args_t>(args)...})::traits_type,
                          typename decltype(sam_file_input{std::forward<args_t>(args)...})::selected_field_ids,
                          typename decltype(sam_file_input{std::forward<args_t>(args)...})::valid_formats>;
//!\endcond
#endif // defined(__cpp_deduction_guides) && __cpp_deduction_guides >= 201907L

} // namespace seqan3

SEQAN3_DEPRECATED_HEADER("This header is deprecated and will be removed in SeqAn-3.1.0 Please #include <seqan3/io/sam_file/input.hpp> instead.")
