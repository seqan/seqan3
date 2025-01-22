// SPDX-FileCopyrightText: 2006-2025 Knut Reinert & Freie Universität Berlin
// SPDX-FileCopyrightText: 2016-2025 Knut Reinert & MPI für molekulare Genetik
// SPDX-License-Identifier: BSD-3-Clause

/*!\file
 * \brief Adds seqan2 to the test environment.
 * \author Marcel Ehrhardt <marcel.ehrhardt AT fu-berlin.de>
 */

#pragma once

#include <seqan3/core/platform.hpp>

// SeqAn 2 [optional]
/*!\def SEQAN3_HAS_SEQAN2
 * \brief Whether SeqAn2 library is available or not.
 * \ingroup test
 */
#if __has_include(<seqan/version.h>)

#    include <seqan/version.h>

// Require minimum version of seqan2, see https://github.com/seqan/seqan3/issues/748
#    if SEQAN_VERSION_MAJOR == 2 && SEQAN_VERSION_MINOR >= 4
#        define SEQAN3_HAS_SEQAN2 1
#    endif

#endif

#ifdef SEQAN3_HAS_SEQAN2

#    include <ranges>

#    include <seqan/basic.h>

namespace std
{

/*!\brief This makes seqan2 iterators ranges-aware, i.e. makes them fulfil std::indirectly_readable which is in turn
 *        needed for std::input_iterator and std::ranges::input_range.
 */
template <typename... args_t>
struct indirectly_readable_traits<seqan2::Iter<args_t...>>
{
    using value_type = typename seqan2::Value<seqan2::Iter<args_t...>>::Type;
};

} // namespace std

#endif // SEQAN3_HAS_SEQAN2
