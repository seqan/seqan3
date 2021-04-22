// -----------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2020, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2020, Knut Reinert & MPI für molekulare Genetik
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at: https://github.com/seqan/seqan3/blob/master/LICENSE.md
// -----------------------------------------------------------------------------------------------------

/*!\file
 * \author Hannes Hauswedell <hannes.hauswedell AT fu-berlin.de>
 * \brief [DEPRECATED] Provides seqan3::views::persist.
 */

#pragma once

#include <seqan3/core/detail/persist_view.hpp>

SEQAN3_DEPRECATED_HEADER(
  "This header is deprecated and will be removed along all its content in SeqAn-3.1.0.")

namespace seqan3::views
{

/*!\brief A view adaptor that wraps rvalue references of non-views.
 * \ingroup views
 * \deprecated This view is deprecated and will be removed in SeqAn3.1.
 */
SEQAN3_DEPRECATED_310 inline auto constexpr persist = detail::persist_fn{};

} // namespace seqan3::views
