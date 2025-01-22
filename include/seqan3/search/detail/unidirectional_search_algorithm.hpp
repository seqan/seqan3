// SPDX-FileCopyrightText: 2006-2025 Knut Reinert & Freie Universität Berlin
// SPDX-FileCopyrightText: 2016-2025 Knut Reinert & MPI für molekulare Genetik
// SPDX-License-Identifier: BSD-3-Clause

/*!\file
 * \author Christopher Pockrandt <christopher.pockrandt AT fu-berlin.de>
 * \brief Provides an approximate string matching algorithm based on simple backtracking.
 *        This should only be used as a reference for unit testing.
 */

#pragma once

#include <ranges>
#include <type_traits>

#include <seqan3/alphabet/concept.hpp>
#include <seqan3/core/detail/test_accessor.hpp>
#include <seqan3/search/detail/search_common.hpp>
#include <seqan3/search/detail/search_traits.hpp>
#include <seqan3/search/fm_index/concept.hpp>

namespace seqan3::detail
{

//!\brief An enumerator for the different error types used during the backtracking.
//!\ingroup search
enum class error_type : uint8_t
{
    deletion,  //!< A deletion was enumerated in previous backtracking step.
    insertion, //!< A insertion was enumerated in previous backtracking step.
    matchmm,   //!< A match or a mismatch was enumerated.
    none       //!< No error or match was enumerated yet.
};

/*!\brief The algorithm that performs a unidirectional search on an FM index using trivial backtracking.
 * \ingroup search
 * \tparam configuration_t The search configuration type.
 * \tparam index_t The type of index.
 */
template <typename configuration_t, typename index_t, typename... policies_t>
class unidirectional_search_algorithm : protected policies_t...
{
private:
    //!\brief The search configuration traits.
    using traits_t = search_traits<configuration_t>;
    //!\brief The search result type.
    using search_result_type = typename traits_t::search_result_type;

    static_assert(!std::same_as<search_result_type, empty_type>, "The search result type was not configured.");

public:
    /*!\name Constructors, destructor and assignment
     * \{
     */
    unidirectional_search_algorithm() = default;                                                    //!< Defaulted.
    unidirectional_search_algorithm(unidirectional_search_algorithm const &) = default;             //!< Defaulted.
    unidirectional_search_algorithm(unidirectional_search_algorithm &&) = default;                  //!< Defaulted.
    unidirectional_search_algorithm & operator=(unidirectional_search_algorithm const &) = default; //!< Defaulted.
    unidirectional_search_algorithm & operator=(unidirectional_search_algorithm &&) = default;      //!< Defaulted.
    ~unidirectional_search_algorithm() = default;                                                   //!< Defaulted.

    /*!\brief Constructs from a configuration object and an index.
     * \tparam configuration_t The search configuration type.
     * \tparam index_t The type of index.
     * \param[in] cfg The configuration object that guides the search algorithm.
     * \param[in] index The index used in the algorithm.
     *
     * \details
     *
     * Initialises the stratum value from the configuration if it was set by the user.
     */
    unidirectional_search_algorithm(configuration_t const & cfg, index_t const & index) : policies_t{cfg}...
    {
        stratum = cfg.get_or(search_cfg::hit_strata{0}).stratum;
        index_ptr = &index;
    }
    //!\}

    /*!\brief Searches a query sequence in an FM index using trivial backtracking.
     *
     * \tparam indexed_query_t The type of the indexed query sequence; must model seqan3::tuple_like with exactly two
     *                         elements and the second tuple element must model std::ranges::forward_range over the
     *                         index's alphabet.
     * \tparam callback_t The callback type to be invoked on a search result; must model std::invocable with the
     *                    search result.
     *
     * \param[in] indexed_query The indexed query sequence to be searched in the index.
     * \param[in] callback The callback to call on a search result.
     *
     * \details
     *
     * The indexed_query parameter is a pair of an index and a query which shall be searched in the index.
     * The search result can then be identified by the index that was associated with the given query.
     *
     * ### Complexity
     *
     * \f$O(|query|^e)\f$ where \f$e\f$ is the maximum number of errors.
     */
    template <typename indexed_query_t, typename callback_t>
        requires (std::tuple_size_v<indexed_query_t> == 2)
              && std::ranges::forward_range<std::tuple_element_t<1, indexed_query_t>>
              && std::invocable<callback_t, search_result_type>
    void operator()(indexed_query_t && indexed_query, callback_t && callback)
    {
        auto && [query_idx, query] = indexed_query;
        auto error_state = this->max_error_counts(query); // see policy_max_error

        // construct internal delegate for collecting hits for later filtering (if necessary)
        std::vector<typename index_t::cursor_type> internal_hits{};
        delegate = [&internal_hits](auto const & it)
        {
            internal_hits.push_back(it);
        };

        perform_search_by_hit_strategy(internal_hits, query, error_state);

        this->make_results(std::move(internal_hits), query_idx, callback); // see policy_search_result_builder
    }

private:
    //!\brief A pointer to the fm index which is used to perform the unidirectional search.
    index_t const * index_ptr{nullptr};

    //!\brief A function object that stores the on-hit-delegate to be executed whenever a hit in the index is found.
    std::function<void(typename index_t::cursor_type const &)> delegate;

    //!\brief The stratum value if set.
    uint8_t stratum{};

    // forward declaration
    template <bool abort_on_hit, typename query_t>
    bool search_trivial(typename index_t::cursor_type cur,
                        query_t & query,
                        typename index_t::cursor_type::size_type const query_pos,
                        search_param const error_left,
                        error_type const prev_error);

    /*!\brief Calls search_trivial depending on the search strategy (hit configuration) given in the configuration.
     * \tparam query_t Must model std::ranges::input_range over the index's alphabet.
     * \param[in, out] internal_hits The result vector to be filled.
     * \param[in] query Query sequence to be searched with the cursor.
     * \param[in] error_state Number of errors for matching the query sequence.
     */
    template <typename query_t>
    void perform_search_by_hit_strategy(std::vector<typename index_t::cursor_type> & internal_hits,
                                        query_t & query,
                                        search_param error_state)
    {
        if constexpr (!traits_t::search_all_hits)
        {
            auto max_total = error_state.total;
            error_state.total = 0; // start search with less errors

            while (internal_hits.empty() && error_state.total <= max_total)
            {
                // * If you only want the best hit (traits_t::search_single_best_hit), you stop after finding the
                //   first hit, the hit with the least errors (`abort_on_hit` is true).
                // * If you are in strata mode (traits_t::search_strata_hits), you do the same as with best hits,
                //   but then do the extra step afterwards (`abort_on_hit` is true).
                // * If you want all best hits (traits_t::search_all_best_hits), you do not stop after the first
                //   hit but continue the current search algorithm/max_error pattern (`abort_on_hit` is true).
                constexpr bool abort_on_hit = !traits_t::search_all_best_hits;
                search_trivial<abort_on_hit>(index_ptr->cursor(), query, 0, error_state, error_type::none);
                ++error_state.total;
            }

            if constexpr (traits_t::search_strata_hits)
            {
                if (!internal_hits.empty())
                {
                    internal_hits.clear();
                    error_state.total += stratum - 1;
                    search_trivial<false>(index_ptr->cursor(), query, 0, error_state, error_type::none);
                }
            }
        }
        else // traits_t::search_all
        {
            // If you want to find all hits, you cannot stop once you found any hit (<false>)
            // since you have to find all paths in the search tree that satisfy the hit condition.
            search_trivial<false>(index_ptr->cursor(), query, 0, error_state, error_type::none);
        }
    }
};

/*!\brief Searches a query sequence in an index using trivial backtracking.
 * \ingroup search
 * \tparam abort_on_hit  If the flag is set, the search algorithm aborts on the first hit.
 * \tparam query_t       Must model std::ranges::input_range over the index's alphabet.
 * \param[in] cur        Cursor of a string index built on the text that will be searched.
 * \param[in] query      Query sequence to be searched with the cursor.
 * \param[in] query_pos  Position in the query sequence indicating the prefix that has already been searched.
 * \param[in] error_left Number of errors left for matching the remaining suffix of the query sequence.
 * \param[in] prev_error Previous scenario of search, i.e. error or match.
 * \returns `True` if and only if `abort_on_hit` is `true` and a hit has been found.
 *
 * ### Complexity
 *
 * \f$O(|query|^e)\f$ where \f$e\f$ is the maximum number of errors.
 *
 * ### Exceptions
 *
 * No-throw guarantee if invoking the delegate also guarantees no-throw.
 */
template <typename configuration_t, typename index_t, typename... policies_t>
template <bool abort_on_hit, typename query_t>
inline bool unidirectional_search_algorithm<configuration_t, index_t, policies_t...>::search_trivial(
    typename index_t::cursor_type cur,
    query_t & query,
    typename index_t::cursor_type::size_type const query_pos,
    search_param const error_left,
    error_type const prev_error)
{
    // Exact case (end of query sequence or no errors left)
    if (query_pos == std::ranges::size(query) || error_left.total == 0)
    {
        // If not at end of query sequence, try searching the remaining suffix without any errors.
        using drop_size_t = std::ranges::range_difference_t<query_t>;
        if (query_pos == std::ranges::size(query)
            || cur.extend_right(std::views::drop(query, static_cast<drop_size_t>(query_pos))))
        {
            delegate(cur);
            return true;
        }
    }
    // Approximate case
    else
    {
        // Insertion
        // Only allow insertions if there is no match and we are not at the beginning of the query.
        bool const allow_insertion =
            (cur.query_length() > 0) ? cur.last_rank() != seqan3::to_rank(query[query_pos]) : true;

        if (allow_insertion && (prev_error != error_type::deletion || error_left.substitution == 0)
            && error_left.insertion > 0)
        {
            search_param error_left2{error_left};
            error_left2.insertion--;
            error_left2.total--;

            // Always perform a recursive call. Abort recursion if and only if recursive call found a hit and
            // abort_on_hit is set to true.
            if (search_trivial<abort_on_hit>(cur, query, query_pos + 1, error_left2, error_type::insertion)
                && abort_on_hit)
            {
                return true;
            }
        }

        // Do not allow deletions at the beginning of the query sequence
        if (((query_pos > 0 && error_left.deletion > 0) || error_left.substitution > 0) && cur.extend_right())
        {
            do
            {
                // Match (when error_left.substitution > 0) and Mismatch
                if (error_left.substitution > 0)
                {
                    bool delta = cur.last_rank() != seqan3::to_rank(query[query_pos]);
                    search_param error_left2{error_left};
                    error_left2.total -= delta;
                    error_left2.substitution -= delta;

                    if (search_trivial<abort_on_hit>(cur, query, query_pos + 1, error_left2, error_type::matchmm)
                        && abort_on_hit)
                    {
                        return true;
                    }
                }

                // Deletion (Do not allow deletions at the beginning of the query sequence.)
                if (query_pos > 0)
                {
                    // Match (when error_left.substitution == 0)
                    if (error_left.substitution == 0 && cur.last_rank() == seqan3::to_rank(query[query_pos]))
                    {
                        if (search_trivial<abort_on_hit>(cur, query, query_pos + 1, error_left, error_type::matchmm)
                            && abort_on_hit)
                        {
                            return true;
                        }
                    }

                    // Deletions at the end of the sequence are not allowed. When the algorithm
                    // arrives here, it cannot be at the end of the query and since deletions do not touch the query
                    // (i.e. increase query_pos) it won't be at the end of the query after the deletion.
                    // Do not allow deletions after an insertion.
                    if ((prev_error != error_type::insertion || error_left.substitution == 0)
                        && error_left.deletion > 0)
                    {
                        search_param error_left2{error_left};
                        error_left2.total--;
                        error_left2.deletion--;
                        // Only search for characters different from the corresponding query character.
                        // (Same character is covered by a match.)
                        if (cur.last_rank() != seqan3::to_rank(query[query_pos]))
                        {
                            if (search_trivial<abort_on_hit>(cur, query, query_pos, error_left2, error_type::deletion)
                                && abort_on_hit)
                            {
                                return true;
                            }
                        }
                    }
                }
            }
            while (cur.cycle_back());
        }
        else
        {
            // Match (when error_left.substitution == 0)
            if (cur.extend_right(query[query_pos]))
            {
                if (search_trivial<abort_on_hit>(cur, query, query_pos + 1, error_left, error_type::matchmm)
                    && abort_on_hit)
                {
                    return true;
                }
            }
        }
    }

    return false;
}

} // namespace seqan3::detail
