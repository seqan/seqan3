// -----------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2020, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2020, Knut Reinert & MPI für molekulare Genetik
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at: https://github.com/seqan/seqan3/blob/master/LICENSE.md
// -----------------------------------------------------------------------------------------------------

/*!\file
 * \author Christopher Pockrandt <christopher.pockrandt AT fu-berlin.de>
 * \brief Provides the public interface for search algorithms.
 */

#pragma once

#include <seqan3/std/algorithm>
#include <seqan3/std/ranges>

#include <seqan3/core/algorithm/configuration.hpp>
#include <seqan3/range/views/persist.hpp>
#include <seqan3/range/views/type_reduce.hpp>
#include <seqan3/search/configuration/default_configuration.hpp>
#include <seqan3/search/detail/policy_max_error.hpp>
#include <seqan3/search/detail/policy_result_builder.hpp>
#include <seqan3/search/detail/unidirectional_search_algorithm.hpp>
#include <seqan3/search/detail/search_scheme_algorithm.hpp>
#include <seqan3/search/detail/search_traits.hpp>
#include <seqan3/search/search_result_range.hpp>

namespace seqan3::detail
{

/*!\brief Class used to update the search configuration, e.g. add defaults.
 * \ingroup search
 */
struct search_configurator
{
    /*!\brief Add seqan3::search_cfg::all to the configuration if seqan3::search_cfg::mode was not set.
     * \tparam configuration_t The type of the search configuration.
     * \param[in] cfg The configuration to be modified if necessary.
     * \returns The configuration which is guaranteed to have a seqan3::search_cfg::mode available.
     *
     * \details
     *
     * If seqan3::search_cfg::mode was not set, it defaults to seqan3::search_cfg::all.
     */
    template <typename configuration_t>
    static auto add_default_mode_configuration(configuration_t const & cfg)
    {
        if constexpr (!detail::search_traits<configuration_t>::has_mode_configuration)
            return cfg | search_cfg::mode{search_cfg::all};
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
     * \sa seqan3::details::search_configurator::add_default_mode_configuration
     * \sa seqan3::details::search_configurator::add_default_output_configuration
     */
    template <typename configuration_t>
    static auto add_defaults(configuration_t const & cfg)
    {
        static_assert(detail::is_type_specialisation_of_v<configuration_t, configuration>,
                      "cfg must be a specialisation of seqan3::configuration.");

        auto cfg1 = add_default_mode_configuration(cfg);
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

/*!\brief Class used to validate the search configuration.
 * \ingroup search
 */
struct search_configuration_validator
{
    /*!\brief Validates the error configuration.
     *
     * \tparam configuration_t The type of the search configuration.
     *
     * \param[in] cfg The configuration to validate.
     *
     * \throws std::invalid_argument
     *
     * \details
     *
     * Checks if the given thresholds for the configured errors are valid. Otherwise throws std::invalid_argument.
     */
    template <typename configuration_t>
    static void validate_error_configuration(configuration_t const & cfg)
    {
        static_assert(detail::is_type_specialisation_of_v<configuration_t, configuration>,
                      "cfg must be a specialisation of seqan3::configuration.");

        using search_traits_t = detail::search_traits<configuration_t>;

        if constexpr (search_traits_t::search_with_max_error)
        {
            auto const & [total, subs, ins, del] = get<search_cfg::max_error>(cfg).value;
            if (subs > total)
                throw std::invalid_argument{"The substitution error threshold is higher than the total error "
                                            "threshold."};
            if (ins > total)
                throw std::invalid_argument{"The insertion error threshold is higher than the total error threshold."};
            if (del > total)
                throw std::invalid_argument{"The deletion error threshold is higher than the total error threshold."};
        }
        else if constexpr (search_traits_t::search_with_max_error_rate)
        {
            auto const & [total, subs, ins, del] = get<search_cfg::max_error_rate>(cfg).value;
            if (subs > total)
                throw std::invalid_argument{"The substitution error threshold is higher than the total error "
                                            "threshold."};
            if (ins > total)
                throw std::invalid_argument{"The insertion error threshold is higher than the total error threshold."};
            if (del > total)
                throw std::invalid_argument{"The deletion error threshold is higher than the total error threshold."};
        }
    }

    /*!\brief Validates the query type to model std::ranges::random_access_range and std::ranges::sized_range.
     *
     * \tparam query_t The type of the query or range of queries.
     */
    template <typename query_t>
    static void validate_query_type()
    {
        using pure_query_t = remove_cvref_t<query_t>;
        if constexpr(dimension_v<pure_query_t> == 1u)
        {
            static_assert(std::ranges::random_access_range<pure_query_t>,
                          "The query sequence must model random_access_range.");
            static_assert(std::ranges::sized_range<pure_query_t>, "The query sequence must model sized_range.");
        }
        else
        {
            static_assert(std::ranges::forward_range<pure_query_t>, "The query collection must model forward_range.");
            static_assert(std::ranges::sized_range<pure_query_t>, "The query collection must model sized_range.");
            static_assert(std::ranges::random_access_range<std::ranges::range_value_t<pure_query_t>>,
                          "Elements of the query collection must model random_access_range.");
            static_assert(std::ranges::sized_range<std::ranges::range_value_t<pure_query_t>>,
                          "Elements of the query collection must model sized_range.");
        }
    }
};
} // namespace seqan3::detail

namespace seqan3
{

/*!\addtogroup search
 * \{
 */

/*!\brief Search a query or a range of queries in an index.
 * \tparam index_t    Must model seqan3::fm_index_specialisation.
 * \tparam queries_t  Must model std::ranges::random_access_range over the index's alphabet and std::ranges::sized_range.
 *                    A range of queries must additionally model std::ranges::forward_range and std::ranges::sized_range.
 * \param[in] queries A single query or a range of queries.
 * \param[in] index   String index to be searched.
 * \param[in] cfg     A configuration object specifying the search parameters (e.g. number of errors, error types,
 *                    output format, etc.).
 * \returns
 * <table>
 *   <tr>
 *     <th>seqan3::text_layout</th>
 *     <th>seqan3::search_cfg::output</th>
 *     <th>Result of seqan3::search()</th>
 *   </tr>
 *   <tr>
 *     <td style="text-align:center">\ref seqan3::text_layout "single"</td>
 *     <td style="text-align:center">\ref seqan3::search_cfg::text_position "text_position"</td>
 *     <td>A `std::vector<size_t>` representing text positions where the search was successful.</td>
 *   </tr>
 *   <tr>
 *     <td style="text-align:center">\ref seqan3::text_layout "single"</td>
 *     <td style="text-align:center">\ref seqan3::search_cfg::index_cursor "index_cursor"</td>
 *     <td>A `std::vector<typename index_t::cursor_type>` containing index_cursors at the text positions where the
 *         search was successful.</td>
 *   </tr>
 *   <tr>
 *     <td style="text-align:center">\ref seqan3::text_layout "collection"</td>
 *     <td style="text-align:center">\ref seqan3::search_cfg::text_position "text_position"</td>
 *     <td>A `std::vector<std::pair<size_t, size_t>>` where the first element of the
 *         `std::pair` specifies the text index in the collection and the second element contains the position in that
 *         text of the collection.</td>
 *   </tr>
 *   <tr>
 *     <td style="text-align:center">\ref seqan3::text_layout "collection"</td>
 *     <td style="text-align:center">\ref seqan3::search_cfg::index_cursor "index_cursor"</td>
 *     <td>A `std::vector<typename index_t::cursor_type>` containing index_cursors at the text positions where the
 *         search was successful.</td>
 *   </tr>
 * </table>
 *
 * \if DEV \note Always returns `void` if an on_hit delegate has been specified.\endif
 *
 * \details
 *
 * \header_file{seqan3/search/search.hpp}
 *
 * ### Complexity
 *
 * Each query with \f$e\f$ errors takes \f$O(|query|^e)\f$ where \f$e\f$ is the maximum number of errors.
 *
 * ### Exceptions
 *
 * Strong exception guarantee if iterating the query does not change its state and if invoking a possible delegate
 * specified in `cfg` also has a strong exception guarantee; basic exception guarantee otherwise.
 *
 * ### Example
 *
 * \include test/snippet/search/search.cpp
 */
template <fm_index_specialisation index_t,
          std::ranges::forward_range queries_t,
          typename configuration_t = decltype(search_cfg::default_configuration)>
//!\cond
    requires std::ranges::forward_range<std::ranges::range_reference_t<queries_t>>
//!\endcond
inline auto search(queries_t && queries,
                   index_t const & index,
                   configuration_t const & cfg = search_cfg::default_configuration)
{
    auto updated_cfg = detail::search_configurator::add_defaults(cfg);

    detail::search_configuration_validator::validate_query_type<queries_t>();
    detail::search_configuration_validator::validate_error_configuration(updated_cfg);

    // return type: for each query: a vector of text_positions (or cursors)
    // delegate params: text_position (or cursor). we will withhold all hits of one query anyway to filter
    //                  duplicates. more efficient to call delegate once with one vector instead of calling
    //                  delegate for each hit separately at once.
    auto algorithm = detail::search_configurator::configure_algorithm(updated_cfg, index);

    return search_result_range{std::move(algorithm), std::forward<queries_t>(queries) | views::type_reduce};
}

//!\cond DEV
// Overload for a single query (not a collection of queries)
//!\overload
template <fm_index_specialisation index_t,
          std::ranges::forward_range query_t,
          typename configuration_t = decltype(search_cfg::default_configuration)>
inline auto search(query_t query,
                   index_t const & index,
                   configuration_t const & cfg = search_cfg::default_configuration)
{
    return search(std::views::single(std::move(query)), index, cfg);
}

//!\overload
template <fm_index_specialisation index_t, typename configuration_t = decltype(search_cfg::default_configuration)>
inline auto search(char const * const queries,
                   index_t const & index,
                   configuration_t const & cfg = search_cfg::default_configuration)
{
    return search(std::string_view{queries}, index, cfg);
}

//!\overload
template <fm_index_specialisation index_t, typename configuration_t = decltype(search_cfg::default_configuration)>
inline auto search(std::initializer_list<char const * const> const & queries,
                   index_t const & index,
                   configuration_t const & cfg = search_cfg::default_configuration)
{
    std::vector<std::string_view> query;
    query.reserve(std::ranges::size(queries));
    std::ranges::for_each(queries, [&query] (char const * const q) { query.push_back(std::string_view{q}); });
    return search(std::move(query) | views::persist, index, cfg);
}
//!\endcond

//!\}

} // namespace seqan3
