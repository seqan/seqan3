// -----------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2018, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2018, Knut Reinert & MPI für molekulare Genetik
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at: https://github.com/seqan/seqan3/blob/master/LICENSE
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

/*!\brief A configuration element for alignment bands.
 * \ingroup configuration
 *
 * \tparam band_t The underlying band class; must be a band data structure.
 */
template <typename band_t>
//!\cond
    requires std::Same<band_t, static_band>
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
