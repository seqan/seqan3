// -----------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2020, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2020, Knut Reinert & MPI für molekulare Genetik
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at: https://github.com/seqan/seqan3/blob/master/LICENSE.md
// -----------------------------------------------------------------------------------------------------

/*!\file
 * \brief [DEPRECATED] Provides seqan3::lower_bound and seqan3::upper_bound.
 * \author Rene Rahn <rene.rahn AT fu-berlin.de>
 * \deprecated This header will be removed in 3.1.0; The contained functionality has been replaced by the
 *             seqan3::align_cfg::band_fixed_size configuration.
 */

#pragma once

#include <seqan3/core/detail/strong_type.hpp>
#include <seqan3/utility/detail/exposition_only_concept.hpp>

SEQAN3_DEPRECATED_HEADER("This header is deprecated and will be removed in SeqAn-3.1.0. The contained functionality has been replaced by the seqan3::align_cfg::band_fixed_size configuration.")

namespace seqan3
{

/*!\brief Type for a lower boundary.
 * \ingroup alignment_band
 * \tparam value_t The underlying type of the lower bound; must model seqan3::arithmetic.
 */
template <seqan3::arithmetic value_t>
struct SEQAN3_DEPRECATED_310 lower_bound : detail::strong_type<value_t, lower_bound<value_t>>
{
    //!\brief Inheriting constructors from base class.
    using detail::strong_type<value_t, lower_bound<value_t>>::strong_type;
};

/*!\brief Type for an upper boundary.
 * \ingroup alignment_band
 * \tparam value_t The underlying type of the upper bound; must model seqan3::arithmetic.
 */
template <seqan3::arithmetic value_t>
struct SEQAN3_DEPRECATED_310 upper_bound : detail::strong_type<value_t, upper_bound<value_t>>
{
    //!\brief Inheriting constructors from base class.
    using detail::strong_type<value_t, upper_bound<value_t>>::strong_type;
};

/*!\name Deduction guides
 * \brief Deduces template parameter from the argument.
 * \{
 */
/*!\brief Deduces the underlying lower boundary type.
 * \relates seqan3::lower_bound
 * \tparam value_t The underlying type of the lower bound; must model seqan3::arithmetic.
 */
template <seqan3::arithmetic value_t>
lower_bound(value_t) -> lower_bound<value_t>;

/*!\brief Deduces the underlying upper boundary type.
 * \tparam value_t The underlying type of the upper bound; must model seqan3::arithmetic.
 */
template <seqan3::arithmetic value_t>
upper_bound(value_t) -> upper_bound<value_t>;
//!\}

} // namespace seqan3
