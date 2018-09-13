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
 * \brief Provides seqan3::detail::align_config_band.
 * \author Jörg Winkler <j.winkler AT fu-berlin.de>
 */

#pragma once

#include <seqan3/alignment/configuration/align_config_band_static.hpp>
#include <seqan3/alignment/configuration/utility.hpp>
#include <seqan3/core/algorithm/all.hpp>
#include <seqan3/core/metafunction/basic.hpp>
#include <seqan3/core/metafunction/template_inspection.hpp>
#include <seqan3/std/concepts>

namespace seqan3::detail
{

/*!\brief Tests whether the passed type can be used in a band configuration element. Defaults to `std::false_type`.
 * \ingroup alignment
 * \see seqan3::detail::is_band_config_v
 *
 * \tparam object_t The type which is tested for being a band configurator.
 */
template <typename object_t>
struct is_band_config : public std::false_type
{};

/*!\brief Helper variable template for seqan3::detail::is_band_config.
 * \ingroup alignment
 * \sa seqan3::detail::is_band_config
 *
 * \tparam object_t The type which is tested for being a band configurator.
 */
template <typename object_t>
inline constexpr bool is_band_config_v = is_band_config<object_t>::value;

/*!\brief A configuration element for alignment bands.
 * \ingroup configuration
 *
 * \tparam band_t The underlying band class; must model seqan3::is_band_config.
 */
template <typename band_t>
//!\cond
    requires is_band_config_v<band_t>
//!\endcond
struct align_config_band
{
    //!\brief Holds the actual band.
    band_t value;
};

/*!\brief The band adaptor enabling pipe notation.
 * \ingroup configuration
 *
 * \tparam band_t A template-template class specifying the actual band implementation.
 */
template <template <typename ...> typename band_t>
//!\cond
    requires is_band_config_v<band_t<int32_t>>
//!\endcond
struct align_config_band_adaptor : public configuration_fn_base<align_config_band_adaptor<band_t>>
{
    /*!\brief Adds to the configuration a band configuration element.
     * \tparam configuration_t The type of the configuration to be extended.
     * \tparam value_t The type for the lower and upper band boundaries; must model std::Integral.
     * \param[in] cfg  The configuration to be extended.
     * \param[in] lower The lower boundary of the band for the banded algorithm.
     * \param[in] upper The upper boundary of the band for the banded algorithm.
     * \returns A new configuration containing the band configuration element.
     */
    template <typename configuration_t, std::Integral value_t>
    //!\cond
        requires is_algorithm_configuration_v<remove_cvref_t<configuration_t>>
    //!\endcond
    constexpr auto invoke(configuration_t && cfg,
                          lower_bound<value_t> const lower,
                          upper_bound<value_t> const upper) const
    {
        static_assert(is_valid_alignment_configuration_v<align_cfg::id::band, remove_cvref_t<configuration_t>>,
                      SEQAN3_INVALID_CONFIG(align_cfg::id::band));

        return std::forward<configuration_t>(cfg).push_front(align_config_band<band_t<value_t>>{{lower, upper}});
    }
};

//!\brief Helper template meta-function associated with seqan3::detail::align_config_band.
//!\ingroup configuration
template <>
struct on_align_config<align_cfg::id::band>
{
    //!\brief Type alias used by meta::find_if
    template <config_element_concept t>
    using invoke = typename is_type_specialisation_of<t, align_config_band>::type;
};

/*!
 * \brief Mapping from the seqan3::detail::align_config_band type to it's corresponding seqan3::align_cfg::id.
 * \ingroup configuration
 *
 * \tparam band_t The underlying band class.
 */
template <typename band_t>
struct align_config_type_to_id<align_config_band<band_t>>
{
    //!\brief The associated seqan3::align_cfg::id.
    static constexpr align_cfg::id value = align_cfg::id::band;
};

//!\cond
template <std::Integral value_t>
struct is_band_config<band_static<value_t>> : public std::true_type
{};
//!\endcond

} // namespace seqan3::detail

namespace seqan3::align_cfg
{
/*!\brief A configuration adaptor for a static band.
 * \ingroup configuration
 */
inline constexpr detail::align_config_band_adaptor<seqan3::band_static> band_static;

} // namespace seqan3::align_cfg
