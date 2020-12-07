// -----------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2020, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2020, Knut Reinert & MPI für molekulare Genetik
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at: https://github.com/seqan/seqan3/blob/master/LICENSE.md
// -----------------------------------------------------------------------------------------------------

/*!\file
 * \brief Provides utility functions for bit twiddling.
 * \author Marcel Ehrhardt <marcel.ehrhardt AT fu-berlin.de>
 * \sa https://en.wikipedia.org/wiki/Bit_manipulation
 */

#pragma once

#include <climits>

#include <seqan3/utility/detail/integer_traits.hpp>

namespace seqan3::detail
{

/*!\brief How many bits has a type?
 * \ingroup utility
 *
 * \tparam type_t The type to determine the number of bits.
 */
template <typename type_t>
constexpr auto bits_of = min_viable_uint_v<CHAR_BIT * sizeof(type_t)>;

} // namespace seqan3::detail
