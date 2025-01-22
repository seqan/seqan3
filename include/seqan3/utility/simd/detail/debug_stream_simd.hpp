// SPDX-FileCopyrightText: 2006-2025 Knut Reinert & Freie Universität Berlin
// SPDX-FileCopyrightText: 2016-2025 Knut Reinert & MPI für molekulare Genetik
// SPDX-License-Identifier: BSD-3-Clause

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

/*!\brief Prints simd types.
 * \tparam simd_t The simd type to print.
 * \ingroup utility_simd
 */
template <simd::simd_concept simd_t>
struct simd_printer<simd_t>
{
    /*!\brief Prints the simd value.
     * \tparam stream_t The type of the stream.
     * \tparam arg_t The type of the argument.
     * \param[in,out] stream The output stream.
     * \param[in] arg The simd type to print.
     */
    template <typename stream_t, typename arg_t>
    constexpr void operator()(stream_t & stream, arg_t && arg) const
    {
        constexpr size_t length = simd::simd_traits<simd_t>::length;
        using scalar_type = typename simd::simd_traits<simd_t>::scalar_type;

        std::array<scalar_type, length> array{};
        for (size_t i = 0; i < length; ++i)
            array[i] = arg[i];
        stream << array;
    }
};

} // namespace seqan3
