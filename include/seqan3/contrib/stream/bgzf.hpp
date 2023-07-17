// -----------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2023, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2023, Knut Reinert & MPI für molekulare Genetik
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at: https://github.com/seqan/seqan3/blob/master/LICENSE.md
// -----------------------------------------------------------------------------------------------------

/*!\file
 * \brief Provides seqan3::contrib::bgzf_thread_count.
 * \author Enrico Seiler <enrico.seiler AT fu-berlin.de>
 */

#pragma once

#include <cstdint>

namespace seqan3::contrib
{

/*!\brief A static variable indicating the number of threads to use for the bgzf-streams. Defaults to 4.
 */
[[maybe_unused]] inline uint64_t bgzf_thread_count = 4;

} // namespace seqan3::contrib
