// -----------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2020, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2020, Knut Reinert & MPI für molekulare Genetik
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at: https://github.com/seqan/seqan3/blob/master/LICENSE.md
// -----------------------------------------------------------------------------------------------------

/*!\file
 * \brief Provides seqan3::detail::align_config_band.
 * \author Jörg Winkler <j.winkler AT fu-berlin.de>
 * \author Rene Rahn <rene.rahn AT fu-berlin.de>
 */

#pragma once

#include <seqan3/alignment/configuration/detail.hpp>
#include <seqan3/alignment/band/static_band.hpp>
#include <seqan3/core/algorithm/pipeable_config_element.hpp>

namespace seqan3::align_cfg
{

/*!\brief Configuration element for setting the band.
 * \ingroup alignment_configuration
 *
 * \tparam band_t The type of the band.
 *
 * \details
 *
 * Configures the banded alignment algorithm. Currently only seqan3::static_band is allowed as argument.
 * If no band is configured for the alignment algorithm the full alignment matrix will be computed.
 * Before executing the algorithm the band is tested for valid settings, e.g. that the upper bound is not smaller than
 * the lower bound, or the band is not shifted out of the alignment matrix. If an invalid setting is detected, a
 * seqan3::invalid_alignment_configuration exception will be thrown.
 *
 * ### Example
 *
 * \include test/snippet/alignment/configuration/align_cfg_band_example.cpp
 */
template <typename band_t>
//!\cond
    requires std::same_as<band_t, static_band>
//!\endcond
struct band : public pipeable_config_element<band<band_t>, band_t>
{
    //!\privatesection
    //!\brief Internal id to check for consistent configuration settings.
    static constexpr detail::align_config_id id{detail::align_config_id::band};
};

/*!\name Type deduction guides
 * \brief Deduces the template parameter from the argument.
 * \relates seqan3::align_cfg::band
 * \{
 */
/*!
 * \brief Deduces the underlying band type.
 * \tparam band_t The underlying type of the band.
 */
template <typename band_t>
band(band_t) -> band<band_t>;
//!\}

} // namespace seqan3::detail
