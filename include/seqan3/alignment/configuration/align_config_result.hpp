// -----------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2018, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2018, Knut Reinert & MPI für molekulare Genetik
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at: https://github.com/seqan/seqan3/blob/master/LICENSE
// -----------------------------------------------------------------------------------------------------

/*!\file
 * \brief Provides configuration for alignment output.
 * \author Rene Rahn <rene.rahn AT fu-berlin.de>
 */

#pragma once

#include <seqan3/alignment/configuration/detail.hpp>
#include <seqan3/alignment/pairwise/align_result.hpp>
#include <seqan3/core/algorithm/pipeable_config_element.hpp>

namespace seqan3::align_cfg
{
/*!\brief A configuration element for specifying the result of the alignment.
 * \ingroup configuration
 */
template <align_result_key e>
class result : public pipeable_config_element
{
public:
    //!\privatesection
    //!\brief An internal id used to check for a valid alignment configuration inherited from the alignment policy.
    static constexpr detail::align_config_id id{detail::align_config_id::result};

    //!\publicsection
    /*!\name Constructor, destructor and assignment
     * \brief Defaulted all standard constructor.
     * \{
     */
    constexpr result()                           noexcept = default;
    constexpr result(result const &)             noexcept = default;
    constexpr result(result &&)                  noexcept = default;
    constexpr result & operator=(result const &) noexcept = default;
    constexpr result & operator=(result &&)      noexcept = default;
    ~result()                                    noexcept = default;
    //!}

    //!\brief The value of align_config_output.
    align_result_key value{e};
};

} //namespace seqan3::align_cfg
