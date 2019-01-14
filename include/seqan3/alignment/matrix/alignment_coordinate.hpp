// -----------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2019, Knut Reinert & Freie Universität Berlin 
// Copyright (c) 2016-2019, Knut Reinert & MPI für molekulare Genetik 
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at: https://github.com/seqan/seqan3/blob/master/LICENSE
// -----------------------------------------------------------------------------------------------------

/*!\file
 * \author Marcel Ehrhardt <marcel.ehrhardt AT fu-berlin.de>
 * \brief Provides seqan3::detail::alignment_coordinate.
 */

#pragma once

#include <seqan3/core/platform.hpp>

namespace seqan3::detail
{

/*!\brief Represents the begin/end of the pairwise alignment in the respective sequences.
 * This class can for example be used to represent the coordinate where the best alignment score is located.
 */
struct alignment_coordinate
{
    //!\brief The position in the first sequence.
    size_t seq1_pos;
    //!\brief The position in the second sequence.
    size_t seq2_pos;
};

} // namespace seqan3::detail
