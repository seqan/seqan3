// -----------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2020, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2020, Knut Reinert & MPI für molekulare Genetik
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at: https://github.com/seqan/seqan3/blob/master/LICENSE.md
// -----------------------------------------------------------------------------------------------------

/*!\file
 * \author Christopher Pockrandt <christopher.pockrandt AT fu-berlin.de>
 * \brief Provides an approximate string matching algorithm based on simple backtracking.
 *        This should only be used as a reference for unit testing.
 */

#pragma once

#include <seqan3/std/ranges>
#include <type_traits>

#include <seqan3/alphabet/concept.hpp>
#include <seqan3/core/detail/test_accessor.hpp>
#include <seqan3/search/detail/search_common.hpp>
#include <seqan3/search/detail/search_traits.hpp>
#include <seqan3/search/fm_index/concept.hpp>
#include <seqan3/range/concept.hpp>
#include <seqan3/range/views/drop.hpp>

namespace seqan3::detail
{

//!\brief An enumerator for the different error types used during the backtracking.
enum class error_type : uint8_t
{
    deletion,  //!< A deletion was enumerated in previous backtracking step.
    insertion, //!< A insertion was enumerated in previous backtracking step.
    matchmm,   //!< A match or a mismatch was enumerated.
    none       //!< No error or match was enumerated yet.
};

/*!\addtogroup search
 * \{
 */

/*!\brief The algorithm that performs a unidirectional search on an FM index using trivial backtracking.
 * \tparam configuration_t The search configuration type.
 * \tparam index_t The type of index; must model seqan3::fm_index_specialisation.
 */
template <typename configuration_t, fm_index_specialisation index_t, typename ...policies_t>
class unidirectional_search_algorithm : protected policies_t...
{
private:
    //!\brief The search configuration traits.
    using traits_t = search_traits<configuration_t>;

public:
    /*!\name Constructors, destructor and assignment
     * \{
     */
    unidirectional_search_algorithm() = default; //!< Defaulted.
    unidirectional_search_algorithm(unidirectional_search_algorithm const &) = default; //!< Defaulted.
    unidirectional_search_algorithm(unidirectional_search_algorithm &&) = default; //!< Defaulted.
    unidirectional_search_algorithm & operator=(unidirectional_search_algorithm const &) = default; //!< Defaulted.
    unidirectional_search_algorithm & operator=(unidirectional_search_algorithm &&) = default; //!< Defaulted.
    ~unidirectional_search_algorithm() = default; //!< Defaulted.

    /*!\brief Constructs from a configuration object and an index.
     * \tparam configuration_t The search configuration type.
     * \tparam index_t The type of index; must model seqan3::fm_index_specialisation.
     * \param[in] cfg The configuration object that guides the search algorithm.
     * \param[in] index The index used in the algorithm.
     */
    unidirectional_search_algorithm(configuration_t const & cfg, index_t const & index) : policies_t{}...
    {
        config = cfg;
        index_ptr = &index;
    }
    //!\}

    /*!\brief Searches a query sequence in an FM index using trivial backtracking.
     * \tparam query_t The type of the query sequence to search; must model std::ranges::input_range over the index's alphabet.
     * \param[in] query Query sequence to be searched in the index.
     *
     * ### Complexity
     *
     * \f$O(|query|^e)\f$ where \f$e\f$ is the maximum number of errors.
     */
    template <typename query_t>
    auto operator()(query_t && query) noexcept
    {
        auto error_state = this->max_error_counts(config, query); // see policy_max_error

        // construct internal delegate for collecting hits for later filtering (if necessary)
        std::vector<typename index_t::cursor_type> internal_hits{};
        delegate = [&internal_hits] (auto const & it)
        {
            internal_hits.push_back(it);
        };

        perform_search_by_mode(internal_hits, query, error_state);

        return this->make_results(std::move(internal_hits), config); // see policy_result_builder
    }

private:
    //!\brief The configuration object.
    configuration_t config{};

    //!\brief A pointer to the fm index which is used to perform the unidirectional search.
    index_t const * index_ptr{nullptr};

    //!\brief A function object that stores the on-hit-delegate to be executed whenever a hit in the index is found.
    std::function<void(typename index_t::cursor_type const &)> delegate;

    // forward declaration
    template <bool abort_on_hit, typename query_t>
    bool search_trivial(typename index_t::cursor_type cur,
                        query_t & query,
                        typename index_t::cursor_type::size_type const query_pos,
                        search_param const error_left,
                        error_type const prev_error);

    /*!\brief Calls search_trivial depending on the search mode given in the configuration.
     * \tparam query_t Must model std::ranges::input_range over the index's alphabet.
     * \param[in, out] internal_hits The result vector to be filled.
     * \param[in] query Query sequence to be searched with the cursor.
     * \param[in] error_state Number of errors for matching the query sequence.
     */
    template <typename query_t>
    void perform_search_by_mode(std::vector<typename index_t::cursor_type> & internal_hits,
                                query_t & query,
                                search_param error_state)
    {
        if constexpr (traits_t::search_best_hits || traits_t::search_all_best_hits || traits_t::search_strata_hits)
        {
            auto max_total = error_state.total;
            error_state.total = 0; // start search with less errors

            while (internal_hits.empty() && error_state.total <= max_total)
            {
                // * If you only want the best hit (traits_t::search_best_hits), you stop after finding the
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
                    uint8_t const stratum = get<search_cfg::mode>(config).value;
                    error_state.total += stratum - 1;
                    search_trivial<false>(index_ptr->cursor(), query, 0, error_state, error_type::none);
                }
            }
        }
        else // detail::search_mode_all
        {
            // If you want to find all hits, you cannot stop once you found any hit (<false>)
            // since you have to find all paths in the search tree that satisfy the hit condition.
            search_trivial<false>(index_ptr->cursor(), query, 0, error_state, error_type::none);
        }
    }

    //!\brief Befriend seqan3::detail::test_accessor to grant access to layout.
    friend struct ::seqan3::detail::test_accessor;
};

/*!\brief Searches a query sequence in an index using trivial backtracking.
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
template <typename configuration_t, typename index_t, typename ...policies_t>
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
        if (query_pos == std::ranges::size(query) || cur.extend_right(views::drop(query, query_pos)))
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
        bool const allow_insertion = (cur.query_length() > 0) ? cur.last_rank() != seqan3::to_rank(query[query_pos]) : true;

        if (allow_insertion && (prev_error != error_type::deletion || error_left.substitution == 0) &&
            error_left.insertion > 0)
        {
            search_param error_left2{error_left};
            error_left2.insertion--;
            error_left2.total--;

            // Always perform a recursive call. Abort recursion if and only if recursive call found a hit and
            // abort_on_hit is set to true.
            if (search_trivial<abort_on_hit>(cur,
                                             query, query_pos + 1,
                                             error_left2,
                                             error_type::insertion) &&
                abort_on_hit)
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

                    if (search_trivial<abort_on_hit>(cur,
                                                     query,
                                                     query_pos + 1,
                                                     error_left2,
                                                     error_type::matchmm) &&
                        abort_on_hit)
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
                        if (search_trivial<abort_on_hit>(cur,
                                                         query,
                                                         query_pos + 1,
                                                         error_left,
                                                         error_type::matchmm) &&
                            abort_on_hit)
                        {
                            return true;
                        }
                    }

                    // Deletions at the end of the sequence are not allowed. When the algorithm
                    // arrives here, it cannot be at the end of the query and since deletions do not touch the query
                    // (i.e. increase query_pos) it won't be at the end of the query after the deletion.
                    // Do not allow deletions after an insertion.
                    if ((prev_error != error_type::insertion || error_left.substitution == 0) &&
                        error_left.deletion > 0)
                    {
                        search_param error_left2{error_left};
                        error_left2.total--;
                        error_left2.deletion--;
                        // Only search for characters different from the corresponding query character.
                        // (Same character is covered by a match.)
                        if (cur.last_rank() != seqan3::to_rank(query[query_pos]))
                        {
                            if (search_trivial<abort_on_hit>(cur,
                                                             query,
                                                             query_pos,
                                                             error_left2,
                                                             error_type::deletion) &&
                                abort_on_hit)
                            {
                                return true;
                            }
                        }
                    }
                }
            } while (cur.cycle_back());
        }
        else
        {
            // Match (when error_left.substitution == 0)
            if (cur.extend_right(query[query_pos]))
            {
                if (search_trivial<abort_on_hit>(cur,
                                                 query,
                                                 query_pos + 1,
                                                 error_left,
                                                 error_type::matchmm) &&
                    abort_on_hit)
                {
                    return true;
                }
            }
        }
    }

    return false;
}
//!\}

} // namespace seqan3::detail
