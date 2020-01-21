// -----------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2020, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2020, Knut Reinert & MPI für molekulare Genetik
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at: https://github.com/seqan/seqan3/blob/master/LICENSE.md
// -----------------------------------------------------------------------------------------------------

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

#include <seqan/version.h>

// Require minimum version of seqan2, see https://github.com/seqan/seqan3/issues/748
#if SEQAN_VERSION_MAJOR == 2 && SEQAN_VERSION_MINOR >= 4
#define SEQAN3_HAS_SEQAN2 1
#endif

#endif

#ifdef SEQAN3_HAS_SEQAN2

#include <seqan/basic.h>
#include <seqan3/std/ranges>

#if __cpp_lib_ranges // C++20 ranges available
namespace std::ranges
#else // implement via range-v3
namespace ranges
#endif
{

//!\brief This makes seqan iterators ranges-aware, i.e. makes them fulfil ranges::readable which is in turn needed for
//!ranges::input_range.
template <typename ...args_t>
struct readable_traits<seqan::Iter<args_t...>>
{
    using value_type = typename seqan::Value<seqan::Iter<args_t...>>::Type;
};

} // namespace ranges
#endif
