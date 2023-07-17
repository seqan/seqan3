// -----------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2023, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2023, Knut Reinert & MPI für molekulare Genetik
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at: https://github.com/seqan/seqan3/blob/master/LICENSE.md
// -----------------------------------------------------------------------------------------------------

/*!\file
 * \author Enrico Seiler <enrico.seiler AT fu-berlin.de>
 * \brief Provides seqan3::views::zip.
 */

#pragma once

#include <seqan3/contrib/std/zip_view.hpp>
#include <seqan3/core/platform.hpp>
#include <seqan3/utility/tuple/common_tuple.hpp> // Included to not break API.

namespace seqan3::views
{

/*!\brief A view adaptor that takes several views and returns tuple-like values from every i-th element of each view.
 * \ingroup utility_views
 * \noapi{This is a implementation of the C++23 zip_view. It will be replaced with std::views::zip.}
 * \sa https://en.cppreference.com/w/cpp/ranges/zip_view
 */
using SEQAN3_DOXYGEN_ONLY(zip =) seqan::std::views::zip;

} // namespace seqan3::views
