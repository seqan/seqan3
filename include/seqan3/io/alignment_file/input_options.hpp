// -----------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2020, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2020, Knut Reinert & MPI für molekulare Genetik
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at: https://github.com/seqan/seqan3/blob/master/LICENSE.md
// -----------------------------------------------------------------------------------------------------

/*!\file
 * \brief [DEPRECATED] Provides seqan3::sam_file_input_options.
 * \deprecated This header will be removed in 3.1.0; Please \#include <seqan3/io/sam_file/input_options.hpp> instead.
 * \author Svenja Mehringer <svenja.mehringer AT fu-berlin.de>
 */

#pragma once

#include <seqan3/io/sam_file/input_options.hpp>

namespace seqan3
{
//!\deprecated Use seqan3::sam_file_input_options instead.
template <typename sequence_legal_alphabet>
using alignment_file_input_options SEQAN3_DEPRECATED_310 = sam_file_input_options<sequence_legal_alphabet>;
} // namespace seqan3

SEQAN3_DEPRECATED_HEADER("This header is deprecated and will be removed in SeqAn-3.1.0 Please #include <seqan3/io/sam_file/input_options.hpp> instead.")
