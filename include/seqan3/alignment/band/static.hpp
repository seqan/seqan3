// ============================================================================
//                 SeqAn - The Library for Sequence Analysis
// ============================================================================
//
// Copyright (c) 2006-2018, Knut Reinert & Freie Universitaet Berlin
// Copyright (c) 2016-2018, Knut Reinert & MPI Molekulare Genetik
// All rights reserved.
//
// Redistribution and use in source and binary forms, with or without
// modification, are permitted provided that the following conditions are met:
//
//     * Redistributions of source code must retain the above copyright
//       notice, this list of conditions and the following disclaimer.
//     * Redistributions in binary form must reproduce the above copyright
//       notice, this list of conditions and the following disclaimer in the
//       documentation and/or other materials provided with the distribution.
//     * Neither the name of Knut Reinert or the FU Berlin nor the names of
//       its contributors may be used to endorse or promote products derived
//       from this software without specific prior written permission.
//
// THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
// AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
// IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
// ARE DISCLAIMED. IN NO EVENT SHALL KNUT REINERT OR THE FU BERLIN BE LIABLE
// FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL
// DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR
// SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER
// CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT
// LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY
// OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH
// DAMAGE.
//
// ============================================================================

/*!\file
 * \brief Provides seqan3::band_static.
 * \author Jörg Winkler <j.winkler AT fu-berlin.de>
 */

#pragma once

#include <seqan3/alignment/band/detail.hpp>
#include <seqan3/core/detail/strong_type.hpp>

namespace seqan3
{

/*!\brief Type for the lower diagonal in a banded alignment.
 *\ingroup alignment
 */
template <typename value_t>
requires std::is_integral_v<value_t>
struct lower_bound : detail::strong_type<value_t, lower_bound<value_t>>
{
    //!\brief Inheriting constructors from base class.
    using detail::strong_type<value_t, lower_bound<value_t>>::strong_type;
};

/*!\brief Type for the upper diagonal in a banded alignment.
 *\ingroup alignment
 */
template <typename value_t>
requires std::is_integral_v<value_t>
struct upper_bound : detail::strong_type<value_t, upper_bound<value_t>>
{
    //!\brief Inheriting constructors from base class.
    using detail::strong_type<value_t, upper_bound<value_t>>::strong_type;
};

/*!\name Deduction guides
 * \brief Deduces template parameter from the argument.
 * \relates seqan3::lower_bound
 * \{
 */
template <typename value_t>
requires std::is_integral_v<value_t>
lower_bound(value_t) -> lower_bound<value_t>;
//!\}


/*!\name Deduction guides
 * \brief Deduces template parameter from the argument.
 * \relates seqan3::upper_bound
 * \{
 */

template <typename value_t>
requires std::is_integral_v<value_t>
upper_bound(value_t) -> upper_bound<value_t>;
//!\}


/*!\brief Data structure for a static band.
 * \ingroup alignment
 * \tparam value_t The value type for the diagonal indices.
 */
template <typename value_t>
requires std::is_integral_v<value_t>
struct band_static
{
    /*!\name Constructor, destructor and assignment
     * \{
     */
    band_static()                                = default;
    band_static(band_static const &)             = default;
    band_static(band_static &&)                  = default;
    band_static & operator=(band_static const &) = default;
    band_static & operator=(band_static &&)      = default;
    ~band_static()                               = default;

    //!\brief Construction from seqan3::lower_bound and seqan3::upper_bound.
    template <typename inner_value_t>
    constexpr band_static(lower_bound<inner_value_t> const lower, upper_bound<inner_value_t> const upper) noexcept
        : lower_bound{std::move(lower.get())}, upper_bound{std::move(upper.get())}
    {}
    //!}

    //!\privatesection
    //!\brief The data member storing the lower and upper diagonal indices.
    value_t lower_bound{std::numeric_limits<value_t>::max()};
    value_t upper_bound{std::numeric_limits<value_t>::max()};
};

/*!\name Deduction guides
 * \brief Deduces template parameter from the argument.
 * \relates seqan3::band_static
 * \{
 */
template <typename value_t>
requires std::is_integral_v<value_t>
band_static(lower_bound<value_t>, upper_bound<value_t>) -> band_static<value_t>;
//!\}

} // namespace seqan3

namespace seqan3::detail
{

//!\cond
template <typename value_t>
requires std::is_integral_v<value_t>
struct is_band_config<band_static<value_t>> : public std::true_type
{};
//!\endcond

} // namespace seqan3::detail
