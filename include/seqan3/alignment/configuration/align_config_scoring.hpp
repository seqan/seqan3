// -----------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2018, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2018, Knut Reinert & MPI für molekulare Genetik
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at: https://github.com/seqan/seqan3/blob/master/LICENSE
// -----------------------------------------------------------------------------------------------------

/*!\file
 * \brief Provides seqan3::align_cfg::scoring.
 * \author Jörg Winkler <j.winkler AT fu-berlin.de>
 */

#pragma once

#include <seqan3/alignment/configuration/detail.hpp>
#include <seqan3/core/metafunction/basic.hpp>
#include <seqan3/alignment/scoring/scoring_scheme_concept.hpp>
#include <seqan3/core/algorithm/pipeable_config_element.hpp>

namespace seqan3::align_cfg
{

/*!\brief The configuration for scoring class.
 * \tparam scoring_scheme_t The type of the scoring scheme. Must satisfy seqan3::scoring_scheme_concept.
 */
template <scoring_scheme_concept<char> scoring_scheme_t>
class scoring : public pipeable_config_element
{
public:
    //!\privatesection
    //!\brief The identifier for this configuration.
    static constexpr detail::align_config_id id{detail::align_config_id::scoring};

    //!\publicsection
    /*!\name Constructor, destructor and assignment
     * \brief Defaulted all standard constructor.
     * \{
     */
    constexpr scoring()                            = default;
    constexpr scoring(scoring const &)             = default;
    constexpr scoring(scoring &&)                  = default;
    constexpr scoring & operator=(scoring const &) = default;
    constexpr scoring & operator=(scoring &&)      = default;
    ~scoring()                                     = default;

    /*!\brief Creates this config from a scoring_scheme.
     * \param scheme The scoring scheme to set. Must satisfy seqan3::scoring_scheme_concept.
     */
    constexpr scoring(scoring_scheme_t scheme) : value{std::move(scheme)}
    {}
    //!}

    //!\brief The stored configuration value.
    scoring_scheme_t value;
};

/*!\names Type deduction guides
 * \relates seqan3::align_cfg::scoring
 * \{
 */

//!\brief Deduces the scoring scheme type from the constructor argument.
template <scoring_scheme_concept<char> scheme_t>
scoring(scheme_t) -> scoring<remove_cvref_t<scheme_t>>;
//!\}

} // namespace seqan3::align_cfg
