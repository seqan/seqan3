// -----------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2020, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2020, Knut Reinert & MPI für molekulare Genetik
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at: https://github.com/seqan/seqan3/blob/master/LICENSE.md
// -----------------------------------------------------------------------------------------------------

/*!\file
 * \brief [DEPRECATED] Provides seqan3::band_static.
 * \author Jörg Winkler <j.winkler AT fu-berlin.de>
 * \deprecated This header will be removed in 3.1.0; The contained functionality has been replaced by the
 *             seqan3::align_cfg::band_fixed_size configuration.
 */

#pragma once

#include <limits>
#include <stdexcept>

#include <seqan3/core/algorithm/bound.hpp>
#include <seqan3/std/concepts>

namespace seqan3
{

/*!\brief [DEPRECATED] Data structure for a static band.
 * \ingroup alignment_band
 * \deprecated use seqan3::align_cfg::band_fixed_size instead.
 */
class SEQAN3_DEPRECATED_310 static_band
{
public:
    constexpr static_band()                                noexcept = default; //!< Defaulted
    constexpr static_band(static_band const &)             noexcept = default; //!< Defaulted
    constexpr static_band(static_band &&)                  noexcept = default; //!< Defaulted
    constexpr static_band & operator=(static_band const &) noexcept = default; //!< Defaulted
    constexpr static_band & operator=(static_band &&)      noexcept = default; //!< Defaulted
    ~static_band()                                         noexcept = default; //!< Defaulted

    /*!\brief Construction from seqan3::lower_bound and seqan3::upper_bound.
     * \tparam input_value_t The input type of the lower and upper band boundaries.
     *
     * \throws std::invalid_argument if upper < lower.
     *
     * \details
     * The boundaries denote the maximum allowed inbalance of insertions and deletions in the alignment.
     * For a symmetric band, choose lower = -upper. The upper boundary must not be smaller than the lower boundary.
     */
    template <std::integral input_value_t>
    constexpr static_band(lower_bound<input_value_t> const lower, upper_bound<input_value_t> const upper)
        : lower_bound{lower.get()}, upper_bound{upper.get()}
    {
        if (lower.get() > upper.get())
        {
            throw std::invalid_argument("The upper boundary must not be smaller than the lower boundary.");
        }
    }

    //!\brief The data member storing the lower boundary of the band.
    int64_t lower_bound{std::numeric_limits<int64_t>::lowest()};
    //!\brief The data member storing the upper boundary of the band.
    int64_t upper_bound{std::numeric_limits<int64_t>::max()};
};

} // namespace seqan3
