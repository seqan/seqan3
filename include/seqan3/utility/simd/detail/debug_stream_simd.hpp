// -----------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2023, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2023, Knut Reinert & MPI für molekulare Genetik
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at: https://github.com/seqan/seqan3/blob/master/LICENSE.md
// -----------------------------------------------------------------------------------------------------

/*!\file
 * \brief Provides seqan3::debug_stream overload for seqan3::simd::simd_type.
 * \author Marcel Ehrhardt <marcel.ehrhardt AT fu-berlin.de>
 */

#pragma once

#include <seqan3/core/debug_stream/debug_stream_type.hpp>
#include <seqan3/core/debug_stream/range.hpp>
#include <seqan3/utility/simd/concept.hpp>
#include <seqan3/utility/simd/simd_traits.hpp>

namespace seqan3
{

/*!\brief Overload for debug_stream for simd types.
 * \ingroup utility_simd
 */
template <typename char_t, typename simd_t>
    requires simd::simd_concept<std::remove_cvref_t<simd_t>>
inline debug_stream_type<char_t> & operator<<(debug_stream_type<char_t> & s, simd_t && simd)
{
    using simd_type = std::remove_cvref_t<simd_t>;
    constexpr size_t length = simd::simd_traits<simd_type>::length;
    using scalar_type = typename simd::simd_traits<simd_type>::scalar_type;

    std::array<scalar_type, length> array{};
    for (size_t i = 0; i < length; ++i)
        array[i] = simd[i];
    s << array;
    return s;
}

} // namespace seqan3
