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
#include <seqan3/core/algorithm/detail/algorithm_executor_blocking.hpp>
#include <seqan3/range/views/persist.hpp>
#include <seqan3/range/views/zip.hpp>
#include <seqan3/search/configuration/default_configuration.hpp>
#include <seqan3/search/detail/search_configurator.hpp>
#include <seqan3/search/detail/search_traits.hpp>
#include <seqan3/search/search_result_range.hpp>

namespace seqan3::detail
{
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
    size_t queries_size = std::ranges::distance(queries);
    auto indexed_queries = views::zip(std::views::iota(size_t{0}, queries_size), queries);

    using indexed_queries_t = decltype(indexed_queries);

    using query_t = std::ranges::range_reference_t<indexed_queries_t>;
    auto [algorithm, complete_config] = detail::search_configurator::configure_algorithm<query_t>(updated_cfg, index);

    using complete_configuration_t = decltype(complete_config);
    using algorithm_result_t = typename detail::search_traits<complete_configuration_t>::search_result_type;
    using execution_handler_t = std::conditional_t<complete_configuration_t::template exists<search_cfg::parallel>(),
                                                   detail::execution_handler_parallel,
                                                   detail::execution_handler_sequential>;
    using executor_t = detail::algorithm_executor_blocking<indexed_queries_t,
                                                           decltype(algorithm),
                                                           algorithm_result_t,
                                                           execution_handler_t>;

    size_t thread_count = complete_config.template value_or<search_cfg::parallel>(0);
    return search_result_range{executor_t{std::move(indexed_queries),
                                          std::move(algorithm),
                                          algorithm_result_t{},
                                          thread_count}};
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
