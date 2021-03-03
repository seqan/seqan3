// -----------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2020, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2020, Knut Reinert & MPI für molekulare Genetik
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at: https://github.com/seqan/seqan3/blob/master/LICENSE.md
// -----------------------------------------------------------------------------------------------------

/*!\file
 * \brief [DEPRECATED] Provides the seqan3::alignment_file_header class.
 * \deprecated This header will be removed in 3.1.0; Please \#include <seqan3/io/sam_file/header.hpp> instead.
 * \author Svenja Mehringer <svenja.mehringer AT fu-berlin.de>
 */

#pragma once

#include <seqan3/io/sam_file/header.hpp>

namespace seqan3
{
//!\deprecated Use seqan3::sam_file_header instead.
template <std::ranges::forward_range ref_ids_type = std::deque<std::string>>
#if defined(__cpp_deduction_guides) && __cpp_deduction_guides >= 201907L
using alignment_file_header SEQAN3_DEPRECATED_310 = sam_file_header<ref_ids_type>;
#else // defined(__cpp_deduction_guides) && __cpp_deduction_guides >= 201907L
// This is a workaround for pre-c++20 as CTAD don't work on type aliases.
struct SEQAN3_DEPRECATED_310 alignment_file_header :
    public sam_file_header<ref_ids_type>
{
//!\cond
    using sam_file_header_t = sam_file_header<ref_ids_type>;
    using sam_file_header_t::sam_file_header_t;
    using sam_file_header_t::operator=;
//!\endcond
};

//!\cond
template <typename ...args_t>
alignment_file_header(args_t && ...args) ->
    alignment_file_header<std::remove_reference_t<decltype(sam_file_header{std::forward<args_t>(args)...}.ref_ids())>>;
//!\endcond
#endif // defined(__cpp_deduction_guides) && __cpp_deduction_guides >= 201907L

} // namespace seqan3

SEQAN3_DEPRECATED_HEADER("This header is deprecated and will be removed in SeqAn-3.1.0 Please #include <seqan3/io/sam_file/header.hpp> instead.")
