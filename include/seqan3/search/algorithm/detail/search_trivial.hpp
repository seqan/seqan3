// -----------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2019, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2019, Knut Reinert & MPI für molekulare Genetik
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at: https://github.com/seqan/seqan3/blob/master/LICENSE.md
// -----------------------------------------------------------------------------------------------------

/*!\file
 * \author Christopher Pockrandt <christopher.pockrandt AT fu-berlin.de>
 * \brief Provides an approximate string matching algorithm based on simple backtracking.
 *        This should only be used as a reference for unit testing.
 */

#pragma once

#include <type_traits>

#include <seqan3/range/concept.hpp>
#include <seqan3/range/view/drop.hpp>
#include <seqan3/search/algorithm/detail/search_common.hpp>
#include <seqan3/std/ranges>

namespace seqan3::detail
{

/*!\addtogroup submodule_search_algorithm
 * \{
 */

/*!\brief Searches a query sequence in an index using trivial backtracking.
 * \tparam abort_on_hit  If the flag is set, the search algorithm aborts on the first hit.
 * \tparam cursor_t      Must model seqan3::FmIndexCursor.
 * \tparam query_t       Must model std::ranges::InputRange over the index's alphabet.
 * \tparam delegate_t    Takes `index::cursor_type` as argument.
 * \param[in] cur        Cursor of atring index built on the text that will be searched.
 * \param[in] query      Query sequence to be searched with the cursor.
 * \param[in] query_pos  Position in the query sequence indicating the prefix that has already been searched.
 * \param[in] error_left Number of errors left for matching the remaining suffix of the query sequence.
 * \param[in] delegate   Function that is called on every hit.
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
template <bool abort_on_hit, typename query_t, typename cursor_t, typename delegate_t>
inline bool search_trivial(cursor_t cur, query_t & query, typename cursor_t::size_type const query_pos,
                           search_param const error_left, delegate_t && delegate) noexcept(noexcept(delegate))
{
    // Exact case (end of query sequence or no errors left)
    if (query_pos == std::ranges::size(query) || error_left.total == 0)
    {
        // If not at end of query sequence, try searching the remaining suffix without any errors.
        if (query_pos == std::ranges::size(query) || cur.extend_right(view::drop(query, query_pos)))
        {
            delegate(cur);
            return true;
        }
    }
    // Approximate case
    else
    {
        // Insertion
        if (error_left.insertion > 0)
        {
            search_param error_left2{error_left};
            error_left2.insertion--;
            error_left2.total--;

            // always perform a recursive call. Abort recursion if and only if recursive call found a hit and
            // abort_on_hit is set to true.
            if (search_trivial<abort_on_hit>(cur, query, query_pos + 1, error_left2, delegate) && abort_on_hit)
                return true;
        }

        // Do not allow deletions at the beginning of the query sequence
        if (((query_pos > 0 && error_left.deletion > 0) || error_left.substitution > 0) && cur.extend_right())
        {
            do
            {
                // Match (when error_left.substitution > 0) and Mismatch
                if (error_left.substitution > 0)
                {
                    bool delta = cur.last_rank() != to_rank(query[query_pos]);
                    search_param error_left2{error_left};
                    error_left2.total -= delta;
                    error_left2.substitution -= delta;

                    if (search_trivial<abort_on_hit>(cur, query, query_pos + 1, error_left2, delegate) && abort_on_hit)
                        return true;
                }

                // Deletion (Do not allow deletions at the beginning of the query sequence.)
                if (query_pos > 0)
                {
                    // Match (when error_left.substitution == 0)
                    if (error_left.substitution == 0 && cur.last_rank() == to_rank(query[query_pos]))
                    {
                        if (search_trivial<abort_on_hit>(cur, query, query_pos + 1, error_left, delegate) &&
                            abort_on_hit)
                        {
                            return true;
                        }
                    }

                    // Deletions at the end of the sequence are not allowed. This cannot happen: when the algorithm
                    // arrives here, it cannot be at the end of the query and since deletions do not touch the query
                    // (i.e. increase query_pos) it won't be at the end of the query after the deletion.
                    if (error_left.deletion > 0)
                    {
                        search_param error_left2{error_left};
                        error_left2.total--;
                        error_left2.deletion--;

                        if (search_trivial<abort_on_hit>(cur, query, query_pos, error_left2, delegate) && abort_on_hit)
                            return true;
                    }
                }
            } while (cur.cycle_back());
        }
        else
        {
            // Match (when error_left.substitution == 0)
            if (cur.extend_right(query[query_pos]))
            {
                if (search_trivial<abort_on_hit>(cur, query, query_pos + 1, error_left, delegate) && abort_on_hit)
                    return true;
            }
        }
    }

    return false;
}

/*!\brief Searches a query sequence in an index using trivial backtracking.
 * \tparam abort_on_hit  If the flag is set, the search algorithm aborts on the first hit.
 * \tparam index_t       Must model seqan3::FmIndex.
 * \tparam query_t       Must model std::ranges::InputRange over the index's alphabet.
 * \tparam delegate_t    Takes `index::cursor_type` as argument.
 * \param[in] index      String index built on the text that will be searched.
 * \param[in] query      Query sequence to be searched in the index.
 * \param[in] error_left Number of errors left for matching the remaining suffix of the query sequence.
 * \param[in] delegate   Function that is called on every hit.
 *
 * ### Complexity
 *
 * \f$O(|query|^e)\f$ where \f$e\f$ is the maximum number of errors.
 *
 * ### Exceptions
 *
 * No-throw guarantee if invoking the delegate also guarantees no-throw.
 */
template <bool abort_on_hit, typename index_t, typename query_t, typename delegate_t>
inline void search_trivial(index_t const & index, query_t & query, search_param const error_left,
                           delegate_t && delegate) noexcept(noexcept(delegate))
{
    search_trivial<abort_on_hit>(index.begin(), query, 0, error_left, delegate);
}

//!\}

} // namespace seqan3::detail
