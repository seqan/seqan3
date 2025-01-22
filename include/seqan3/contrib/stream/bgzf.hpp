// SPDX-FileCopyrightText: 2006-2025 Knut Reinert & Freie Universität Berlin
// SPDX-FileCopyrightText: 2016-2025 Knut Reinert & MPI für molekulare Genetik
// SPDX-License-Identifier: BSD-3-Clause

/*!\file
 * \brief Provides seqan3::contrib::bgzf_thread_count.
 * \author Enrico Seiler <enrico.seiler AT fu-berlin.de>
 */

#pragma once

#include <cstdint>

namespace seqan3::contrib
{

/*!\brief A global variable indicating the number of threads to use for the bgzf-streams. Defaults to 4.
 */
[[maybe_unused]] inline uint64_t bgzf_thread_count = 4;

} // namespace seqan3::contrib
