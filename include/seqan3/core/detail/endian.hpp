// -----------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2020, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2020, Knut Reinert & MPI für molekulare Genetik
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at: https://github.com/seqan/seqan3/blob/master/LICENSE.md
// -----------------------------------------------------------------------------------------------------

/*!\file
 * \brief Provides seqan3::detail::endian.
 * \author Rene Rahn <rene.rahn AT fu-berlin.de>
 */

#pragma once

#include <seqan3/core/platform.hpp>

namespace seqan3::detail
{

/*!\brief Indicates the endianness of all scalar types
 * \ingroup core
 *
 * \details
 *
 * This class reimplements the std::endian from https://en.cppreference.com/w/cpp/types/endian available in c++20.
 */
enum class endian
{
#ifdef _WIN32
    little = 0,                       //!< implementation-defined
    big    = 1,                       //!< implementation-defined
    native = little                   //!< implementation-defined
#else
    little = __ORDER_LITTLE_ENDIAN__, //!< implementation-defined
    big    = __ORDER_BIG_ENDIAN__,    //!< implementation-defined
    native = __BYTE_ORDER__           //!< implementation-defined
#endif
};

} // namespace seqan3::detail
