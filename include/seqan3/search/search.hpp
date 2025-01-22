// SPDX-FileCopyrightText: 2006-2025 Knut Reinert & Freie Universität Berlin
// SPDX-FileCopyrightText: 2016-2025 Knut Reinert & MPI für molekulare Genetik
// SPDX-License-Identifier: BSD-3-Clause

/*!\file
 * \brief Provides the public interface for search algorithms.
 * \author Christopher Pockrandt <christopher.pockrandt AT fu-berlin.de>
 */

#pragma once

#include <algorithm>
#include <ranges>

#include <seqan3/core/algorithm/algorithm_result_generator_range.hpp>
#include <seqan3/core/algorithm/detail/algorithm_executor_blocking.hpp>
#include <seqan3/core/configuration/configuration.hpp>
#include <seqan3/search/configuration/default_configuration.hpp>
#include <seqan3/search/configuration/on_result.hpp>
#include <seqan3/search/configuration/parallel.hpp>
#include <seqan3/search/detail/search_configurator.hpp>
#include <seqan3/search/detail/search_traits.hpp>
#include <seqan3/utility/views/convert.hpp>
#include <seqan3/utility/views/deep.hpp>
#include <seqan3/utility/views/zip.hpp>

namespace seqan3::detail
{
/*!\brief Class used to validate the search configuration.
 * \ingroup search
 */
struct search_configuration_validator
{
    /*!\brief Validates the query type to model std::ranges::random_access_range and std::ranges::sized_range.
     *
     * \tparam query_t The type of the query or range of queries.
     */
    template <typename query_t>
    static void validate_query_type()
    {
        using pure_query_t = std::remove_cvref_t<query_t>;
        if constexpr (range_dimension_v<pure_query_t> == 1u)
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
/*!\brief Search a query or a range of queries in an index.
 * \ingroup search
 * \tparam index_t    Type of the index. See \ref search_available_indices for an overview of our indices.
 * \tparam queries_t  Must model std::ranges::random_access_range over the index's alphabet and std::ranges::sized_range.
 *                    A range of queries must additionally model std::ranges::forward_range and std::ranges::sized_range.
 * \param[in] queries A single query or a range of queries.
 * \param[in] index   String index to be searched.
 * \param[in] cfg     A configuration object specifying the search parameters (e.g. number of errors, error types,
 *                    output format, etc.).
 * \returns A seqan3::algorithm_result_generator_range with value type of seqan3::search_result.
 *
 * \if DEV \note Always returns `void` if an on_hit delegate has been specified.\endif
 *
 * \details
 *
 * \header_file{seqan3/search/search.hpp}
 *
 * The search algorithm strongly depends on the **index** that is used.
 * Please see \ref search_available_indices for an overview of our indices.
 *
 * For more details on how to configure the search, please see the respective documentation: \ref search_configuration.
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
template <typename index_t,
          std::ranges::forward_range queries_t,
          typename configuration_t = decltype(search_cfg::default_configuration)>
    requires std::ranges::forward_range<std::ranges::range_reference_t<queries_t>>
          && std::same_as<range_innermost_value_t<queries_t>, typename index_t::alphabet_type>
inline auto
search(queries_t && queries, index_t const & index, configuration_t const & cfg = search_cfg::default_configuration)
{
    auto updated_cfg = detail::search_configurator::add_defaults(cfg);

    detail::search_configuration_validator::validate_query_type<queries_t>();

    size_t queries_size = std::ranges::distance(queries);
    auto indexed_queries = views::zip(std::views::iota(size_t{0}, queries_size), std::forward<queries_t>(queries));

    using indexed_queries_t = decltype(indexed_queries);

    using query_t = std::ranges::range_reference_t<indexed_queries_t>;
    auto [algorithm, complete_config] = detail::search_configurator::configure_algorithm<query_t>(updated_cfg, index);

    using complete_configuration_t = decltype(complete_config);
    using traits_t = detail::search_traits<complete_configuration_t>;
    using algorithm_result_t = typename traits_t::search_result_type;
    using execution_handler_t = std::conditional_t<complete_configuration_t::template exists<search_cfg::parallel>(),
                                                   detail::execution_handler_parallel,
                                                   detail::execution_handler_sequential>;

    // Select the execution handler for the search configuration.
    auto select_execution_handler = [parallel = complete_config.get_or(search_cfg::parallel{})]()
    {
        if constexpr (std::same_as<execution_handler_t, detail::execution_handler_parallel>)
        {
            auto thread_count = parallel.thread_count;
            if (!thread_count)
                throw std::runtime_error{"You must configure the number of threads in seqan3::search_cfg::parallel."};

            return execution_handler_t{*thread_count};
        }
        else
        {
            return execution_handler_t{};
        }
    };

    // Finally, choose between two way execution returning an algorithm range or calling a user callback on every hit.
    if constexpr (traits_t::has_user_callback)
    {
        select_execution_handler().bulk_execute(algorithm,
                                                indexed_queries,
                                                get<search_cfg::on_result>(complete_config).callback);
    }
    else
    {
        using executor_t = detail::algorithm_executor_blocking<indexed_queries_t,
                                                               decltype(algorithm),
                                                               algorithm_result_t,
                                                               execution_handler_t>;

        return algorithm_result_generator_range{executor_t{std::move(indexed_queries),
                                                           std::move(algorithm),
                                                           algorithm_result_t{},
                                                           select_execution_handler()}};
    }
}

//!\cond DEV
// Convert query sequence if it does not match the alphabet type of the index.
//!\overload
template <typename index_t,
          std::ranges::forward_range queries_t,
          typename configuration_t = decltype(search_cfg::default_configuration)>
    requires std::ranges::forward_range<std::ranges::range_reference_t<queries_t>>
          && (!std::same_as<range_innermost_value_t<queries_t>, typename index_t::alphabet_type>)
inline auto search(queries_t && queries,
                   index_t const & index,
                   configuration_t const & cfg = search_cfg::default_configuration)
{
    static_assert(std::convertible_to<range_innermost_value_t<queries_t>, typename index_t::alphabet_type>,
                  "The alphabet of the text collection must be convertible to the alphabet of the index.");

    if constexpr (range_dimension_v<queries_t> == 2u)
        return search(queries | views::deep{views::convert<typename index_t::alphabet_type>}, index, cfg);
    else
        return search(queries | views::convert<typename index_t::alphabet_type>, index, cfg);
}

// Overload for a single query (not a collection of queries)
//!\overload
template <typename index_t,
          std::ranges::forward_range query_t,
          typename configuration_t = decltype(search_cfg::default_configuration)>
inline auto
search(query_t && query, index_t const & index, configuration_t const & cfg = search_cfg::default_configuration)
{
    return search(std::views::single(std::forward<query_t>(query)), index, cfg);
}

//!\overload
template <typename index_t, typename configuration_t = decltype(search_cfg::default_configuration)>
inline auto search(char const * const queries,
                   index_t const & index,
                   configuration_t const & cfg = search_cfg::default_configuration)
{
    return search(std::string_view{queries}, index, cfg);
}

//!\overload
template <typename index_t, typename configuration_t = decltype(search_cfg::default_configuration)>
inline auto search(std::initializer_list<char const * const> const & queries,
                   index_t const & index,
                   configuration_t const & cfg = search_cfg::default_configuration)
{
    std::vector<std::string_view> query;
    query.reserve(std::ranges::size(queries));
    std::ranges::for_each(queries,
                          [&query](char const * const q)
                          {
                              query.push_back(std::string_view{q});
                          });
    return search(std::move(query) | std::views::all, index, cfg);
}
//!\endcond

} // namespace seqan3
