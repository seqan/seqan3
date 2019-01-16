// -----------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2019, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2019, Knut Reinert & MPI für molekulare Genetik
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at: https://github.com/seqan/seqan3/blob/master/LICENSE
// -----------------------------------------------------------------------------------------------------

/*!\file
 * \brief This header includes C++17 filesystem support and imports it into namespace seqan3::filesystem (independent of whether it is marked as "experimental").
 * \author Vitor C. Piro <pirov AT zedat.fu-berlin.de >
 */

#pragma once

#if __has_include(<filesystem>)
#include <filesystem>
#else
#include <experimental/filesystem>
#endif // __has_include(experimental/filesystem)

#include <seqan3/core/platform.hpp>

namespace seqan3
{
#if __has_include(<filesystem>)
namespace filesystem = std::filesystem;
#else
namespace filesystem = std::experimental::filesystem;
#endif // __has_include(experimental/filesystem)
} // namespace seqan3
