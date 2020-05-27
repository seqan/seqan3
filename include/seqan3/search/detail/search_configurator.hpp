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

#include <seqan3/core/type_list/type_list.hpp>
#include <seqan3/core/type_traits/function.hpp>
#include <seqan3/core/type_traits/lazy.hpp>
#include <seqan3/core/type_traits/template_inspection.hpp>
#include <seqan3/search/configuration/hit.hpp>
#include <seqan3/search/configuration/max_error.hpp>
#include <seqan3/search/configuration/result_type.hpp>
#include <seqan3/search/detail/policy_max_error.hpp>
#include <seqan3/search/detail/policy_search_result_builder.hpp>
#include <seqan3/search/detail/search_scheme_algorithm.hpp>
#include <seqan3/search/detail/unidirectional_search_algorithm.hpp>

namespace seqan3::detail
{

/*!\brief Class used to update the search configuration, e.g. add defaults.
 * \ingroup search
 */
class search_configurator
{
private:
    /*!\brief Select the search result based on the configuration and the index type.
     *
     * \tparam search_configuration_t The type of the configuration.
     * \tparam index_t The type of the index.
     */
    template <typename search_configuration_t, typename index_t>
    struct select_search_result
    {
    private:
        //!\brief The position type based on the index text layout.
        using position_type = std::conditional_t<index_t::text_layout_mode == text_layout::collection,
                                                 std::pair<typename index_t::size_type, typename index_t::size_type>,
                                                 typename index_t::size_type>;

    public:
        //!\brief The result type depending on the output configuration.
        using type = std::conditional_t<search_traits<search_configuration_t>::search_return_index_cursor,
                                        typename index_t::cursor_type,
                                        position_type>;
    };

    /*!\brief Selects the search algorithm based on the index type.
     *
     * \tparam search_configuration_t The type of the configuration.
     * \tparam index_t The type of the index.
     * \tparam policies_t A template parameter pack over the policies to specify the behavior of the algorithm.
     */
    template <typename configuration_t, typename index_t, typename ...policies_t>
    struct select_search_algorithm
    {
        //!\brief The selected algorithm type based on the index.
        using type =
            lazy_conditional_t<bi_fm_index_specialisation<index_t>,
                               lazy<instantiate_t,
                                    lazy<search_scheme_algorithm, configuration_t, index_t, policies_t...>>,
                               lazy<instantiate_t,
                                    lazy<unidirectional_search_algorithm, configuration_t, index_t, policies_t...>>>;
    };

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
     *
     * \tparam query_t An explicit template argument for the query type the search algorithm is invoked with.
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
    template <typename query_t, typename configuration_t, typename index_t>
    static auto configure_algorithm(configuration_t const & cfg, index_t const & index)
    {
        using search_result_t = typename select_search_result<configuration_t, index_t>::type;
        using search_result_collection_t = std::vector<search_result_t>;
        using type_erased_algorithm_t = std::function<search_result_collection_t(query_t)>;

        return configure_hit_strategy<type_erased_algorithm_t>(cfg | search_cfg::detail::result_type<search_result_t>,
                                                               index);
    }

    template <typename algorithm_t, typename configuration_t, typename index_t>
    static algorithm_t configure_hit_strategy(configuration_t const &, index_t const &);
};

/*!\brief Configures the algorithm with the correct hit strategy.
 *
 * \tparam algorithm_t The type erased algorithm used for the fixed return type.
 * \tparam configuration_t The type of the search configuration.
 * \tparam index_t The type of the index.
 *
 * \returns The configured search algorithm.
 *
 * \details
 *
 * If the algorithm was configured with the dynamic hit configuration element seqan3::search_cfg::hit, the
 * configuration element is removed and replaced by the selected static hit configuration element.
 * If the hit configuration element is already a static one nothing is changed in the configuration.
 * After selecting the correct hit strategy the corresponding search algorithm is created with the new configuration
 * and the given index.
 *
 * \throws std::invalid_argument if the dynamic hit configuration was not initialised with a hit strategy.
 *
 * If no hit configuration was configured a static assert is emitted during compilation.
 */
template <typename algorithm_t, typename configuration_t, typename index_t>
algorithm_t search_configurator::configure_hit_strategy(configuration_t const & cfg, index_t const & index)
{
    // Select the actual search algorithm and return an instance constructed with the adapted config and the index.
    auto generate_algorithm = [&] (auto new_cfg) -> algorithm_t
    {
        using new_configuration_t = decltype(new_cfg);
        using selected_algorithm_t =
            typename select_search_algorithm<new_configuration_t,
                                             index_t,
                                             policy_max_error,
                                             policy_search_result_builder<new_configuration_t>>::type;

        return selected_algorithm_t{new_cfg, index};
    };

    // Check if dynamic config present, otherwise continue.
    if constexpr (configuration_t::template exists<search_cfg::hit>())
    {
        auto hit_variant = get<search_cfg::hit>(cfg).value;

        if (std::holds_alternative<empty_type>(hit_variant))
            throw std::invalid_argument{"The dynamic hit strategy was not initialised! "
                                        "Please refer to the configuration documentation of the search algorithm for "
                                        "more details."};

        // Remove dynamic config first.
        auto cfg_without_hit = cfg.template remove<search_cfg::hit>();

        // Apply the correct static configuration element.
        return std::visit(multi_invocable
        {
            [&] (hit_all_best_tag) { return generate_algorithm(cfg_without_hit | search_cfg::hit_all_best); },
            [&] (hit_single_best_tag) { return generate_algorithm(cfg_without_hit | search_cfg::hit_single_best); },
            [&] (search_cfg::hit_strata const & strata) { return generate_algorithm(cfg_without_hit | strata); },
            [&] (auto) { return generate_algorithm(cfg_without_hit | search_cfg::hit_all); }
        }, hit_variant);
    }
    else // Already statically configured.
    {
        static_assert(detail::search_traits<configuration_t>::has_hit_configuration,
                      "The hit strategy for the search algorithm was not configured. "
                      "Please refer to the configuration documentation of the search algorithm for more details.");

        return generate_algorithm(cfg);
    }
}

} // namespace seqan3::detail
