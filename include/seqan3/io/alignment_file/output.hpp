// -----------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2020, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2020, Knut Reinert & MPI für molekulare Genetik
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at: https://github.com/seqan/seqan3/blob/master/LICENSE.md
// -----------------------------------------------------------------------------------------------------

/*!\file
 * \brief [DEPRECATED] Provides seqan3::sam_file_output and corresponding traits classes.
 * \deprecated This header will be removed in 3.1.0; Please \#include <seqan3/io/sam_file/output.hpp> instead.
 * \author Svenja Mehringer <svenja.mehringer AT fu-berlin.de>
 */

#pragma once

#include <seqan3/io/sam_file/output.hpp>

#include <seqan3/io/alignment_file/output_options.hpp>
#include <seqan3/io/alignment_file/header.hpp>

namespace seqan3
{
//!\deprecated Use seqan3::sam_file_output instead.
template <detail::fields_specialisation selected_field_ids_ =
              fields<field::seq,
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
          detail::type_list_of_sam_file_output_formats valid_formats_ = type_list<format_sam, format_bam>,
          typename ref_ids_type = ref_info_not_given>
#if defined(__cpp_deduction_guides) && __cpp_deduction_guides >= 201907L
using alignment_file_output SEQAN3_DEPRECATED_310 = sam_file_output<selected_field_ids_, valid_formats_, ref_ids_type>;
#else // defined(__cpp_deduction_guides) && __cpp_deduction_guides >= 201907L
// This is a workaround for pre-c++20 as CTAD don't work on type aliases.
struct SEQAN3_DEPRECATED_310 alignment_file_output :
    public sam_file_output<selected_field_ids_, valid_formats_, ref_ids_type>
{
//!\cond
    using sam_file_output_t = sam_file_output<selected_field_ids_, valid_formats_, ref_ids_type>;
    using sam_file_output_t::sam_file_output_t;
    using sam_file_output_t::operator=;
//!\endcond
};

//!\cond
template <typename ...args_t>
alignment_file_output(args_t && ...args) ->
    alignment_file_output<typename decltype(sam_file_output{std::forward<args_t>(args)...})::selected_field_ids,
                          typename decltype(sam_file_output{std::forward<args_t>(args)...})::valid_formats,
                          typename decltype(sam_file_output{std::forward<args_t>(args)...})::ref_ids_type_t>;
//!\endcond
#endif // defined(__cpp_deduction_guides) && __cpp_deduction_guides >= 201907L

} // namespace seqan3

SEQAN3_DEPRECATED_HEADER("This header is deprecated and will be removed in SeqAn-3.1.0 Please #include <seqan3/io/sam_file/output.hpp> instead.")
