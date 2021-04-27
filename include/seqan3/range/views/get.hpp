// -----------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2020, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2020, Knut Reinert & MPI für molekulare Genetik
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at: https://github.com/seqan/seqan3/blob/master/LICENSE.md
// -----------------------------------------------------------------------------------------------------

/*!\file
 * \author Joshua Kim <joshua.kim AT fu-berlin.de>
 * \brief [DEPRECATED] Provides seqan3::views::get.
 * \deprecated This header will be removed in 3.1. Please use seqan3/utility/views/elements.hpp instead.
 */

#pragma once

#include <seqan3/utility/views/elements.hpp>

namespace seqan3::views
{

/*!\brief A view calling `get` on each element in a range.
 * \ingroup views
 * \deprecated Please use `seqan3::views::elements` instead.
 */
template <auto index>
SEQAN3_DEPRECATED_310 inline constexpr auto get = views::elements<index>;

} // namespace seqan3::views

SEQAN3_DEPRECATED_HEADER(
    "This header is deprecated and will be removed in SeqAn-3.1.0; Please #include <seqan3/utility/views/elements.hpp> instead.")
