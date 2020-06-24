// -----------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2020, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2020, Knut Reinert & MPI für molekulare Genetik
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at: https://github.com/seqan/seqan3/blob/master/LICENSE.md
// -----------------------------------------------------------------------------------------------------

/*!\file
 * \brief Provides global and local alignment configurations.
 * \author Joshua Kim <joshua.kim AT fu-berlin.de>
 * \author Rene Rahn <rene.rahn AT fu-berlin.de>
 * \author Jörg Winkler <j.winkler AT fu-berlin.de>
 * \author Wiep van der Toorn <w.vandertoorn AT fu-berlin.de>
 */

#pragma once

#include <seqan3/alignment/configuration/detail.hpp>
#include <seqan3/core/algorithm/pipeable_config_element.hpp>
#include <seqan3/core/detail/empty_type.hpp>

namespace seqan3::detail
{
    //!\brief A strong type to select the global alignment method.
    //!\ingroup alignment_configuration
    struct method_global_tag : public pipeable_config_element<method_global_tag>
    {
        //!\privatesection
        //!\brief An internal id used to check for a valid alignment configuration.
        static constexpr detail::align_config_id id{detail::align_config_id::global};
    };

    //!\brief A strong type to select the local alignment method.
    //!\ingroup alignment_configuration
    struct method_local_tag : public pipeable_config_element<method_local_tag>
    {
        //!\privatesection
        //!\brief An internal id used to check for a valid alignment configuration.
        static constexpr detail::align_config_id id{detail::align_config_id::local};
    };
} // namespace seqan3::detail

namespace seqan3::align_cfg
{
    /*!\brief Sets the global alignment method.
     * \ingroup alignment_configuration
     *
     * \details
     *
     * The alignment algorithm can be categorised in different methods. For example, the
     * \ref seqan3::align_cfg::method_local "local" and the
     * \ref seqan3::align_cfg::method_global "global" alignment are two different methods, while the semi-global alignment
     * is a variation of the global alignment. This differentiation makes it possible to define a subset of configurations
     * that can work with a particular method. Since it is not possible to guess what the desired method for a user is, this
     * configuration must be provided for the alignment algorithm and cannot be defaulted.
     *
     * ### Example
     *
     * \include test/snippet/alignment/configuration/minimal_alignment_config.cpp
     */
    inline constexpr seqan3::detail::method_global_tag method_global{};

    /*!\brief Sets the local alignment method.
     * \ingroup alignment_configuration
     * \copydetails seqan3::align_cfg::method_global
     */
    inline constexpr seqan3::detail::method_local_tag method_local{};
} // namespace seqan3::align_cfg
