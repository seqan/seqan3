// -----------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2020, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2020, Knut Reinert & MPI für molekulare Genetik
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at: https://github.com/seqan/seqan3/blob/master/LICENSE.md
// -----------------------------------------------------------------------------------------------------

/*!\file
 * \brief [DEPRECATED] Provides seqan3::alignment_file_output_format and auxiliary classes.
 * \deprecated This header will be removed in 3.1.0; Please \#include <seqan3/io/sam_file/output_format_concept.hpp> instead.
 * \author Svenja Mehringer <svenja.mehringer AT fu-berlin.de>
 */

#pragma once

#include <seqan3/io/sam_file/output_format_concept.hpp>

namespace seqan3
{
//!\deprecated Use seqan3::sam_file_output_format instead.
template <typename t>
SEQAN3_DEPRECATED_310 constexpr bool alignment_file_output_format = sam_file_output_format<t>;
} // namespace seqan3

SEQAN3_DEPRECATED_HEADER("This header is deprecated and will be removed in SeqAn-3.1.0 Please #include <seqan3/io/sam_file/output_format_concept.hpp> instead.")
