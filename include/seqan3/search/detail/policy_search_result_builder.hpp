// SPDX-FileCopyrightText: 2006-2025 Knut Reinert & Freie Universität Berlin
// SPDX-FileCopyrightText: 2016-2025 Knut Reinert & MPI für molekulare Genetik
// SPDX-License-Identifier: BSD-3-Clause

/*!\file
 * \author Svenja Mehringer <svenja.mehringer AT fu-berlin.de>
 * \brief Provides the seqan3::detail::policy_search_result_builder.
 */

#pragma once

#include <seqan3/core/detail/template_inspection.hpp>
#include <seqan3/search/detail/search_traits.hpp>
#include <seqan3/search/fm_index/concept.hpp>
#include <seqan3/search/search_result.hpp>

namespace seqan3::detail
{

//!\brief Provides the function `make_results` if inherited by a search algorithm.
//!\ingroup search
template <typename search_configuration_t>
    requires seqan3::detail::is_type_specialisation_of_v<search_configuration_t, configuration>
struct policy_search_result_builder
{
protected:
    //!\brief The traits type over the search configuration.
    using search_traits_type = detail::search_traits<search_configuration_t>;
    //!\brief The configured search result type.
    using search_result_type = typename search_traits_type::search_result_type;

    static_assert(!std::same_as<search_result_type, typename search_traits_type::empty_search_result_type>,
                  "The search result type was not configured properly.");

    /*!\name Constructors, destructor and assignment
     * \{
     */
    policy_search_result_builder() = default;                                                 //!< Defaulted.
    policy_search_result_builder(policy_search_result_builder &&) = default;                  //!< Defaulted.
    policy_search_result_builder(policy_search_result_builder const &) = default;             //!< Defaulted.
    policy_search_result_builder & operator=(policy_search_result_builder &&) = default;      //!< Defaulted.
    policy_search_result_builder & operator=(policy_search_result_builder const &) = default; //!< Defaulted.
    ~policy_search_result_builder() = default;                                                //!< Defaulted.

    //!\brief Construction from the configuration object.
    explicit policy_search_result_builder(search_configuration_t const &)
    {}
    //!\}

    /*!\brief Invoke the callback on all hits (index cursors) without calling locate on each cursor.
     *
     * \tparam index_cursor_t The type of index cursor used in the search algorithm.
     * \tparam query_index_t The index type of the query.
     * \tparam callback_t The callback which is called for every hit.
     *
     * \param[in] internal_hits internal_hits A range over internal cursor results.
     * \param[in] idx The index associated with the current query.
     * \param[in] callback The callback to invoke for every hit.
     *
     * \details
     *
     * The result is independent from the search modus (all, single_best, all_best, strata).
     */
    template <typename index_cursor_t, typename query_index_t, typename callback_t>
    void make_results(std::vector<index_cursor_t> internal_hits, query_index_t idx, callback_t && callback)
    {
        return make_results_impl(std::move(internal_hits), idx, std::forward<callback_t>(callback));
    }

    /*!\brief Invokes the callback on each seqan3::search_result after calling locate on each cursor.
     *
     * \tparam index_cursor_t The type of index cursor used in the search algorithm.
     * \tparam query_index_t The index type of the query.
     * \tparam callback_t The callback which is called for every hit.
     *
     * \param[in] internal_hits internal_hits A range over internal cursor results.
     * \param[in] idx The index associated with the current query.
     * \param[in] callback The callback to invoke for every hit.
     *
     * \details
     *
     * This function is used for all search modi except single_best (which are all, all_best, and strata).
     *
     * The text positions are sorted and made unique by position before invoking the callback on them.
     */
    template <typename index_cursor_t, typename query_index_t, typename callback_t>
        requires search_traits_type::output_requires_locate_call && (!search_traits_type::search_single_best_hit)
    void make_results(std::vector<index_cursor_t> internal_hits, query_index_t idx, callback_t && callback)
    {
        std::vector<search_result_type> results{};
        results.reserve(internal_hits.size()); // expect at least as many text positions as cursors, possibly more

        make_results_impl(std::move(internal_hits),
                          idx,
                          [&results](auto && search_result)
                          {
                              results.push_back(std::move(search_result));
                          });

        // sort by reference id or by reference position if both have the same reference id.
        std::sort(results.begin(),
                  results.end(),
                  [](auto const & r1, auto const & r2)
                  {
                      return (r1.reference_id() == r2.reference_id())
                               ? (r1.reference_begin_position() < r2.reference_begin_position())
                               : (r1.reference_id() < r2.reference_id());
                  });

        results.erase(std::unique(results.begin(), results.end()), results.end());

        for (auto && search_result : results)
            callback(std::move(search_result));
    }

private:
    /*!\brief Invokes the callback on each seqan3::search_result and calls locate on the cursor depending on the config.
     *
     * \tparam index_cursor_t The type of index cursor used in the search algorithm.
     * \tparam query_index_t The index type of the query.
     * \tparam callback_t The callback which is called for every hit.
     *
     * \param[in] internal_hits internal_hits A range over internal cursor results.
     * \param[in] idx The index associated with the current query.
     * \param[in] callback The callback to invoke for every hit.
     *
     * \details
     *
     * For each cursor `in internal_hits`, this function calls `cursor.layz_locate()` if the search configuration
     * requires it (search_traits_type::output_requires_locate_call) and then constructs a seqan3::search_result from
     * the resulting data. The seqan3::search_result will be filled only with the data that was asked for by the user
     * via the `search_traits_type::output_[...]` trait (e.g. `search_traits_type::output_query_id`).
     */
    template <typename index_cursor_t, typename query_index_t, typename callback_t>
    void make_results_impl(std::vector<index_cursor_t> internal_hits,
                           [[maybe_unused]] query_index_t idx,
                           callback_t && callback)
    {
        auto maybe_locate = [](auto const & cursor)
        {
            if constexpr (search_traits_type::output_requires_locate_call)
                return cursor.lazy_locate();
            else
                return std::views::single(std::tuple{0, 0});
        };

        for (auto const & cursor : internal_hits)
        {
            for (auto && [ref_id, ref_pos] : maybe_locate(cursor))
            {
                search_result_type result{};

                if constexpr (search_traits_type::output_query_id)
                    result.query_id_ = std::move(idx);
                if constexpr (search_traits_type::output_index_cursor)
                    result.cursor_ = cursor;
                if constexpr (search_traits_type::output_reference_id)
                    result.reference_id_ = std::move(ref_id);
                if constexpr (search_traits_type::output_reference_begin_position)
                    result.reference_begin_position_ = std::move(ref_pos);

                callback(result);

                if constexpr (search_traits_type::search_single_best_hit)
                    return;
            }
        }
    }
};

} // namespace seqan3::detail
