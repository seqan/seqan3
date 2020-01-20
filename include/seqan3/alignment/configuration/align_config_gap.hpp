// -----------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2020, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2020, Knut Reinert & MPI für molekulare Genetik
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at: https://github.com/seqan/seqan3/blob/master/LICENSE.md
// -----------------------------------------------------------------------------------------------------

/*!\file
 * \brief Provides seqan3::align_config::gap.
 * \author Rene Rahn <rene.rahn AT fu-berlin.de>
 * \author Jörg Winkler <j.winkler AT fu-berlin.de>
 */

#pragma once

#include <seqan3/alignment/configuration/detail.hpp>
#include <seqan3/core/type_traits/basic.hpp>
#include <seqan3/core/type_traits/template_inspection.hpp>
#include <seqan3/alignment/scoring/gap_scheme.hpp>
#include <seqan3/core/algorithm/pipeable_config_element.hpp>

namespace seqan3::align_cfg
{
/*!\brief A configuration element for the gap scheme.
 * \ingroup alignment_configuration
 * \tparam gap_scheme_t The type of the underlying gap scheme; must be of type seqan3::gap_scheme.
 *
 * \details
 *
 * Configures the gap scheme for the alignment algorithm. The gap scheme determines how gaps are penalised inside
 * of the alignment algorithm. If the gap scheme is not configured, it will default to a linear gap scheme initialised
 * with edit distance. Note that the gap open score is used as an additional score. This means that the score for
 * opening a gap during the affine alignment execution is the sum of the gap score and the gap open score.
 *
 * ### Example
 *
 * \include test/snippet/alignment/configuration/align_cfg_gap_example.cpp
 *
 * \see seqan3::gap_scheme
 */
template <typename gap_scheme_t>
struct gap : public pipeable_config_element<gap<gap_scheme_t>, gap_scheme_t>
{
    static_assert(detail::is_type_specialisation_of_v<gap_scheme_t, gap_scheme>,
                  "Expects seqan3::gap_scheme class.");
    //!\privatesection
    //!\brief Internal id to check for consistent configuration settings.
    static constexpr detail::align_config_id id{detail::align_config_id::gap};
};

/*!\name Type deduction guides
 * \relates seqan3::align_cfg::gap
 * \{
 */
//!\brief Deduces the gap scheme from the constructor argument.
template <typename scheme_t>
gap(scheme_t) -> gap<scheme_t>;
//!\}

} // namespace seqan3::align_cfg
