// -----------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2020, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2020, Knut Reinert & MPI für molekulare Genetik
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at: https://github.com/seqan/seqan3/blob/master/LICENSE.md
// -----------------------------------------------------------------------------------------------------

/*!\file
 * \brief Provides seqan3::detail::search_configurator.
 * \author Rene Rahn <rene.rahn AT fu-berlin.de>
 */

#pragma once

#include <seqan3/search/detail/policy_max_error.hpp>
#include <seqan3/search/detail/policy_result_builder.hpp>
#include <seqan3/search/detail/search_scheme_algorithm.hpp>
#include <seqan3/search/detail/unidirectional_search_algorithm.hpp>

namespace seqan3::detail
{

/*!\brief Class used to update the search configuration, e.g. add defaults.
 * \ingroup search
 */
class search_configurator
{
public:
    /*!\brief Add seqan3::search_cfg::hit_all to the configuration if no search strategy (hit configuration) was chosen.
     * \tparam configuration_t The type of the search configuration.
     * \param[in] cfg The configuration to be modified if necessary.
     * \returns The configuration which is guaranteed to have a hit configuration element available.
     *
     * \details
     *
     * If no \ref search_configuration_subsection_hit_strategy "hit configuration" was set,
     * it defaults to seqan3::search_cfg::hit_all.
     */
    template <typename configuration_t>
    static auto add_default_hit_configuration(configuration_t const & cfg)
    {
        if constexpr (!detail::search_traits<configuration_t>::has_hit_configuration)
            return cfg | search_cfg::hit_all;
        else
            return cfg;
    }

    /*!\brief Add seqan3::search_cfg::text_position to the configuration if seqan3::search_cfg::output was not set.
     * \tparam configuration_t The type of the search configuration.
     * \param[in] cfg The configuration to be modified if necessary.
     * \returns The configuration which is guaranteed to have a seqan3::search_cfg::output available.
     *
     * \details
     *
     * If seqan3::search_cfg::output was not set, it defaults to seqan3::search_cfg::text_position.
     */
    template <typename configuration_t>
    static auto add_default_output_configuration(configuration_t const & cfg)
    {
        if constexpr (!detail::search_traits<configuration_t>::has_output_configuration)
            return cfg | search_cfg::output{search_cfg::text_position};
        else
            return cfg;
    }

    /*!\brief Adds default configurations if they were not set by the user.
     * \tparam configuration_t The type of the search configuration.
     * \param[in] cfg The configuration to be modified.
     * \returns The modified configuration.
     *
     * \details
     *
     * Modifies the configuration object by adding default configuration elements.
     *
     * \sa seqan3::details::search_configurator::add_default_hit_configuration
     * \sa seqan3::details::search_configurator::add_default_output_configuration
     */
    template <typename configuration_t>
    static auto add_defaults(configuration_t const & cfg)
    {
        static_assert(detail::is_type_specialisation_of_v<configuration_t, configuration>,
                      "cfg must be a specialisation of seqan3::configuration.");

        auto cfg1 = add_default_hit_configuration(cfg);
        auto cfg2 = add_default_output_configuration(cfg1);

        return cfg2;
    }

    /*!\brief Chooses the appropriate search algorithm depending on the index.
     * \tparam configuration_t The type of the search configuration.
     * \tparam index_t The type of the index.
     * \param[in] cfg The search configuration object that is passed to the algorithm.
     * \param[in] index The index that is passed to the algorithm.
     * \returns A search algorithm.
     *
     * \details
     *
     * If the `index_t` models seqan3::bi_fm_index_specialisation, then the
     * seqan3::detail::search_scheme_algorithm is chosen. Otherwise, the
     * detail::unidirectional_search_algorithm is chosen.
     */
    template <typename configuration_t, typename index_t>
    static auto configure_algorithm(configuration_t const & cfg, index_t const & index)
    {
        if constexpr (bi_fm_index_specialisation<index_t>)
        {
            using algorithm_t = search_scheme_algorithm<configuration_t,
                                                        index_t,
                                                        policy_max_error,
                                                        policy_result_builder>;
            return algorithm_t{cfg, index};
        }
        else
        {
            using algorithm_t = unidirectional_search_algorithm<configuration_t,
                                                                index_t,
                                                                policy_max_error,
                                                                policy_result_builder>;
            return algorithm_t{cfg, index};
        }
    }
};

} // namespace seqan3::detail
