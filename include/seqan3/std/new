// -*- C++ -*-
// -----------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2022, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2022, Knut Reinert & MPI für molekulare Genetik
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at: https://github.com/seqan/seqan3/blob/master/LICENSE.md
// -----------------------------------------------------------------------------------------------------

/*!\file
 * \brief The [\<new\> header](https://en.cppreference.com/w/cpp/header/new) from C++17's standard library.
 * \author Rene Rahn <rene.rahn AT fu-berlin.de>
 */

// File might be included from multiple libraries.
#ifndef SEQAN_STD_NEW_SHIM
#define SEQAN_STD_NEW_SHIM

#include <new>

#ifndef __cpp_lib_hardware_interference_size

/*!\defgroup std_new new
 * \ingroup std
 * \brief The [\<new\> header](https://en.cppreference.com/w/cpp/header/new) from C++17's standard library.
 */

namespace std
{

/*!\brief Minimum offset between two objects to avoid false sharing.
 * \ingroup std_new
 * \sa https://en.cppreference.com/w/cpp/thread/hardware_destructive_interference_size
 */
inline constexpr std::size_t hardware_destructive_interference_size = 64;

/*!\brief Maximum size of contiguous memory to promote true sharing.
 * \ingroup std_new
 * \sa https://en.cppreference.com/w/cpp/thread/hardware_destructive_interference_size
 */
inline constexpr std::size_t hardware_constructive_interference_size = 64;

} // namespace std

#endif // __cpp_lib_hardware_interference_size

#endif // SEQAN_STD_NEW_SHIM
