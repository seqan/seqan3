// -----------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2020, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2020, Knut Reinert & MPI für molekulare Genetik
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at: https://github.com/seqan/seqan3/blob/master/LICENSE.md
// -----------------------------------------------------------------------------------------------------

/*!\file
 * \author Hannes Hauswedell <hannes.hauswedell AT fu-berlin.de>
 * \brief [DEPRECATED] Provides seqan3::views::take_line and seqan3::views::take_line_or_throw.
 */

#pragma once

#include <seqan3/io/detail/take_line_view.hpp>

SEQAN3_DEPRECATED_HEADER(
  "This header is deprecated and will be removed along all its content in SeqAn-3.1.0.")

namespace seqan3::views
{

/*!\brief A view adaptor that returns a single line from the underlying range.
 * \ingroup views
 * \deprecated Please use `std::views::take_while([](auto && chr){ return chr != &apos;\n&apos; })`
 */
SEQAN3_DEPRECATED_310 inline auto constexpr take_line = detail::take_line;

/*!\brief A view adaptor that returns a single line from the underlying range (throws if there is no end-of-line).
 * \throws seqan3::unexpected_end_of_input If the underlying range contains no end-of-line marker.
 * \ingroup views
 * \deprecated Please use `std::views::take_while([](auto && chr){ return chr != &apos;\n&apos; })`
 */
SEQAN3_DEPRECATED_310 inline auto constexpr take_line_or_throw = detail::take_line_or_throw;

} // namespace seqan3::views
