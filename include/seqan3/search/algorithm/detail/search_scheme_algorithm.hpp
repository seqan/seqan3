// -----------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2020, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2020, Knut Reinert & MPI für molekulare Genetik
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at: https://github.com/seqan/seqan3/blob/master/LICENSE.md
// -----------------------------------------------------------------------------------------------------

/*!\file
 * \author Christopher Pockrandt <christopher.pockrandt AT fu-berlin.de>
 * \brief Provides the algorithm to search in an index using search schemes.
 */

#pragma once

#include <type_traits>

#include <seqan3/core/type_traits/transformation_trait_or.hpp>
#include <seqan3/range/views/slice.hpp>
#include <seqan3/search/algorithm/detail/search_common.hpp>
#include <seqan3/search/algorithm/detail/search_scheme_precomputed.hpp>
#include <seqan3/search/algorithm/detail/search_trivial.hpp>
#include <seqan3/search/fm_index/concept.hpp>

namespace seqan3::detail
{

/*!\addtogroup submodule_search_algorithm
 * \{
 */

/*!\brief Computes a (non-optimal) search scheme. Currently the generated search scheme represents trivial backtracking.
 * \param[in] min_error Minimum number of errors allowed.
 * \param[in] max_error Maximum number of errors allowed.
 *
 * ### Complexity
 *
 * Constant.
 *
 * ### Exceptions
 *
 * Strong exception guarantee.
 */
inline std::vector<search_dyn> compute_ss(uint8_t const min_error, uint8_t const max_error)
{
    // TODO: Replace this at least by the pigeonhole principle or even better by 01*0 schemes.
    // NOTE: Make sure that the searches are sorted by their asymptotical running time (i.e. upper error bound string),
    //       s.t. easy to compute searches come first. This improves the running time of algorithms that abort after the
    //       first hit (e.g. search mode: best). Even though it is not guaranteed, this seems to be a good greedy
    //       approach.
    std::vector<search_dyn> scheme{{{1}, {min_error}, {max_error}}};
    return scheme;
}

/*!\brief Returns for each search the cumulative length of blocks in the order of blocks in each search and the
 *        starting position of the first block in the query sequence.
 * \tparam search_scheme_t  Is of type `seqan3::detail::search_scheme_type` or `seqan3::detail::search_scheme_dyn_type`.
 * \param[in] search_scheme Search scheme that will be used for searching.
 * \param[in] query_length  Length of the query that will be searched in an index.
 * \returns A range of pairs containing for each search the cumulative lengths of blocks and the starting position
 *          in the query.
 *
 * ### Complexity
 *
 * Constant.
 *
 * ### Exceptions
 *
 * Strong exception guarantee.
 */
template <typename search_scheme_t>
inline auto search_scheme_block_info(search_scheme_t const & search_scheme, size_t const query_length)
{
    using blocks_length_type = typename search_scheme_t::value_type::blocks_length_type;

    bool constexpr is_dyn_scheme = std::same_as<search_scheme_t, search_scheme_dyn_type>;

    // Either store information in an array (for search schemes known at compile time) or in a vector otherwise.
    using result_type = std::conditional_t<is_dyn_scheme,
                                           std::vector<std::tuple<blocks_length_type, size_t>>,
                                           std::array<std::tuple<blocks_length_type, size_t>,
                                                      transformation_trait_or_t<std::tuple_size<search_scheme_t>,
                                                                                std::false_type>::value>>;

    result_type result;
    if constexpr (is_dyn_scheme)
        result.resize(search_scheme.size());

    uint8_t const blocks      {search_scheme[0].blocks()};
    size_t  const block_length{query_length / blocks};
    uint8_t const rest        {static_cast<uint8_t>(query_length % blocks)};

    blocks_length_type blocks_length;
    // set all blocks_length values to block_length
    // resp. block_length + 1 for the first `rest = block_length % blocks` values
    if constexpr (is_dyn_scheme)
        blocks_length.resize(blocks, block_length);
    else
        blocks_length.fill(block_length);

    for (uint8_t block_id = 0; block_id < rest; ++block_id)
        ++blocks_length[block_id];

    for (uint8_t search_id = 0; search_id < search_scheme.size(); ++search_id)
    {
        auto const & search = search_scheme[search_id];

        auto & [search_blocks_length, start_pos] = result[search_id];

        // compute cumulative blocks_length and starting position
        start_pos = 0;
        if constexpr (is_dyn_scheme)
            search_blocks_length.resize(blocks);
        search_blocks_length[0] = blocks_length[search.pi[0] - 1];
        for (uint8_t i = 1; i < blocks; ++i)
        {
            search_blocks_length[i] = blocks_length[search.pi[i] - 1] + search_blocks_length[i - 1];
            if (search.pi[i] < search.pi[0])
                start_pos += search_blocks_length[i] - search_blocks_length[i - 1];
        }
    }

    return result;
}

//!\cond
// forward declaration
template <bool abort_on_hit, typename cursor_t, typename query_t, typename search_t, typename blocks_length_t,
          typename delegate_t>
inline bool search_ss(cursor_t cur, query_t & query,
                      typename cursor_t::size_type const lb, typename cursor_t::size_type const rb,
                      uint8_t const errors_spent, uint8_t const block_id, bool const go_right, search_t const & search,
                      blocks_length_t const & blocks_length, search_param const error_left, delegate_t && delegate);
//!\endcond

/*!\brief Searches a query sequence in a bidirectional index using a single search of a search scheme.
 *        Sub-function for searching the remaining part of the current block without any errors.
 * \tparam abort_on_hit     If the flag is set, the search aborts on the first hit.
 * \tparam cursor_t         Must model seqan3::bi_fm_index_cursor_specialisation.
 * \tparam query_t          Must model std::ranges::random_access_range over the index's alphabet.
 * \tparam search_t         Is of type `seqan3::detail::search<>` or `seqan3::detail::search_dyn<>`.
 * \tparam blocks_length_t  Is of type `std::array` or `std::vector` of unsigned integers.
 * \tparam delegate_t       Takes `cursor_t` as argument.
 * \param[in] cur           Cursor of a string index built on the text that will be searched.
 * \param[in] query         Query sequence to be searched.
 * \param[in] lb            Left bound of the infix of `query` already searched (exclusive).
 * \param[in] rb            Right bound of the infix of `query` already searched (exclusive).
 * \param[in] errors_spent  Number of errors spent while searching the infix of `query`.
 * \param[in] block_id      Id of the block that the infix is extended to next.
 * \param[in] go_right      The infix will be extended to the right if the flag is set to true.
 * \param[in] search        Search of a search scheme to be used for searching.
 * \param[in] blocks_length Cumulative block lengths of the search.
 * \param[in] error_left    Number of errors left for matching the remaining suffix of the query sequence.
 * \param[in] delegate      Function that is called on every hit.
 * \returns `True` if and only if `abort_on_hit` is true and a hit has been found.
 *
 * ### Complexity
 *
 * \f$O(|query|^e)\f$ where \f$e\f$ is the total number of errors allowed by `search`.
 *
 * ### Exceptions
 *
 * Strong exception guarantee if iterating the query does not change its state and if invoking the delegate also has a
 * strong exception guarantee; basic exception guarantee otherwise.
 */
template <bool abort_on_hit, typename cursor_t, typename query_t, typename search_t, typename blocks_length_t,
          typename delegate_t>
inline bool search_ss_exact(cursor_t cur, query_t & query,
                            typename cursor_t::size_type const lb, typename cursor_t::size_type const rb,
                            uint8_t const errors_spent, uint8_t const block_id, bool const go_right,
                            search_t const & search, blocks_length_t const & blocks_length,
                            search_param const error_left, delegate_t && delegate)
{
    using size_type = typename cursor_t::size_type;

    uint8_t const block_id2 = std::min<uint8_t>(block_id + 1, search.blocks() - 1);
    bool const go_right2 = (block_id < search.blocks() - 1) && (search.pi[block_id + 1] > search.pi[block_id]);

    if (go_right)
    {
        size_type const infix_lb = rb - 1; // inclusive
        size_type const infix_rb = lb + blocks_length[block_id] - 1; // exclusive

        if (!cur.extend_right(query | views::slice(infix_lb, infix_rb + 1)))
            return false;

        if (search_ss<abort_on_hit>(cur, query, lb, infix_rb + 2, errors_spent, block_id2, go_right2, search,
                                    blocks_length, error_left, delegate) && abort_on_hit)
        {
            return true;
        }
    }
    else
    {
        size_type const infix_lb = rb - blocks_length[block_id] - 1; // inclusive
        size_type const infix_rb = lb - 1; // inclusive

        if (!cur.extend_left(query | views::slice(infix_lb, infix_rb + 1)))
            return false;

        if (search_ss<abort_on_hit>(cur, query, infix_lb, rb, errors_spent, block_id2, go_right2, search, blocks_length,
                                    error_left, delegate) && abort_on_hit)
        {
            return true;
        }
    }
    return false;
}

/*!\brief Searches a query sequence in a bidirectional index using a single search of a search schemes.
 *        Sub-function for deletions at the end of a block.
 *
 * \copydetails search_ss_exact
 */
template <bool abort_on_hit, typename cursor_t, typename query_t, typename search_t, typename blocks_length_t,
          typename delegate_t>
inline bool search_ss_deletion(cursor_t cur, query_t & query,
                               typename cursor_t::size_type const lb, typename cursor_t::size_type const rb,
                               uint8_t const errors_spent, uint8_t const block_id, bool const go_right,
                               search_t const & search, blocks_length_t const & blocks_length,
                               search_param const error_left, delegate_t && delegate)
{
    uint8_t const max_error_left_in_block = search.u[block_id] - errors_spent;
    uint8_t const min_error_left_in_block = std::max(search.l[block_id] - errors_spent, 0);

    // Switch to the next block when the min number of errors is reached
    if (min_error_left_in_block == 0)
    {
        uint8_t const block_id2 = std::min<uint8_t>(block_id + 1, search.blocks() - 1);
        bool const go_right2 = search.pi[block_id2] > search.pi[block_id2 - 1];

        if (search_ss<abort_on_hit>(cur, query, lb, rb, errors_spent, block_id2, go_right2, search, blocks_length,
                                    error_left, delegate) && abort_on_hit)
        {
            return true;
        }
    }

    // Insert deletions into the current block as long as possible
    // Do not allow deletions at the beginning of the leftmost block
    // Do not allow deletions at the end of the rightmost block
    if (!(search.pi[block_id] == 1 && !go_right) &&
        !(search.pi[block_id] == search.blocks() && go_right) &&
        max_error_left_in_block > 0 && error_left.total > 0 && error_left.deletion > 0 &&
        ((go_right && cur.extend_right()) || (!go_right && cur.extend_left())))
    {
        search_param error_left2{error_left};
        error_left2.total--;
        error_left2.deletion--;
        do
        {
            if (search_ss_deletion<abort_on_hit>(cur, query, lb, rb, errors_spent + 1, block_id, go_right, search,
                                                 blocks_length, error_left2, delegate) && abort_on_hit)
            {
                return true;
            }
        } while ((go_right && cur.cycle_back()) || (!go_right && cur.cycle_front()));
    }
    return false;
}

/*!\brief Searches a query sequence in a bidirectional index using a single search of a search schemes.
 *        Sub-function for approximate search step (iterating over all children in a conceptual suffix tree).
 *
 * \copydetails search_ss_exact
 *
 * \param[in] min_error_left_in_block Number of remaining errors that need to be spent in the current block.
 */
template <bool abort_on_hit, typename cursor_t, typename query_t, typename search_t, typename blocks_length_t,
          typename delegate_t>
inline bool search_ss_children(cursor_t cur, query_t & query,
                               typename cursor_t::size_type const lb, typename cursor_t::size_type const rb,
                               uint8_t const errors_spent, uint8_t const block_id, bool const go_right,
                               uint8_t const min_error_left_in_block, search_t const & search,
                               blocks_length_t const & blocks_length, search_param const error_left,
                               delegate_t && delegate)
{
    using size_type = typename cursor_t::size_type;
    if ((go_right && cur.extend_right()) || (!go_right && cur.extend_left()))
    {
        size_type const chars_left = blocks_length[block_id] - (rb - lb - 1);

        size_type lb2 = lb - !go_right;
        size_type rb2 = rb + go_right;

        do
        {
            bool const delta = cur.last_rank() != to_rank(query[(go_right ? rb : lb) - 1]);

            // skip if there are more min errors left in the current block than characters in the block
            // i.e. chars_left - 1 < min_error_left_in_block - delta
            // TODO: move that outside the if / do-while struct
            // TODO: incorporate error_left.deletion into formula
            if (error_left.deletion == 0 && chars_left + delta < min_error_left_in_block + 1u)
                continue;

            if (!delta || error_left.substitution > 0)
            {
                search_param error_left2{error_left};
                error_left2.total -= delta;
                error_left2.substitution -= delta;

                // At the end of the current block
                if (rb - lb == blocks_length[block_id])
                {
                    // Leave the possibility for one or multiple deletions at the end of a block.
                    // Thus do not change the direction (go_right) yet.
                    if (error_left.deletion > 0)
                    {
                        if (search_ss_deletion<abort_on_hit>(cur, query, lb2, rb2, errors_spent + delta, block_id,
                                                             go_right, search, blocks_length, error_left2, delegate) &&
                            abort_on_hit)
                        {
                            return true;
                        }
                    }
                    else
                    {
                        uint8_t const block_id2 = std::min<uint8_t>(block_id + 1, search.blocks() - 1);
                        bool const go_right2 = search.pi[block_id2] > search.pi[block_id2 - 1];

                        if (search_ss<abort_on_hit>(cur, query, lb2, rb2, errors_spent + delta, block_id2, go_right2,
                                                    search, blocks_length, error_left2, delegate) &&
                            abort_on_hit)
                        {
                            return true;
                        }
                    }
                }
                else
                {
                    if (search_ss<abort_on_hit>(cur, query, lb2, rb2, errors_spent + delta, block_id, go_right, search,
                                                blocks_length, error_left2, delegate) && abort_on_hit)
                    {
                        return true;
                    }
                }
            }

            // Deletion
            // TODO: check whether the conditions for deletions at the beginning/end of the query are really necessary
            // No deletion at the beginning of the leftmost block.
            // No deletion at the end of the rightmost block.
            if (error_left.deletion > 0 &&
                !(go_right && (rb == 1 || rb == std::ranges::size(query) + 1)) &&
                !(!go_right && (lb == 0 || lb == std::ranges::size(query))))
            {
                search_param error_left3{error_left};
                error_left3.total--;
                error_left3.deletion--;
                search_ss<abort_on_hit>(cur, query, lb, rb, errors_spent + 1, block_id, go_right, search, blocks_length,
                                        error_left3, delegate);
            }
        } while ((go_right && cur.cycle_back()) || (!go_right && cur.cycle_front()));
    }
    return false;
}

/*!\brief Searches a query sequence in a bidirectional index using a single search of a search schemes.
 *
 * \copydetails search_ss_exact
 */
template <bool abort_on_hit, typename cursor_t, typename query_t, typename search_t,
          typename blocks_length_t, typename delegate_t>
inline bool search_ss(cursor_t cur, query_t & query,
                      typename cursor_t::size_type const lb, typename cursor_t::size_type const rb,
                      uint8_t const errors_spent, uint8_t const block_id, bool const go_right, search_t const & search,
                      blocks_length_t const & blocks_length, search_param const error_left, delegate_t && delegate)
{
    uint8_t const max_error_left_in_block = search.u[block_id] - errors_spent;
    uint8_t const min_error_left_in_block = std::max(search.l[block_id] - errors_spent, 0); // NOTE: changed

    // Done.
    if (min_error_left_in_block == 0 && lb == 0 && rb == std::ranges::size(query) + 1)
    {
        delegate(cur);
        return true;
    }
    // Exact search in current block.
    else if (((max_error_left_in_block == 0) && (rb - lb - 1 != blocks_length[block_id])) ||
             (error_left.total == 0 && min_error_left_in_block == 0))
    {
        if (search_ss_exact<abort_on_hit>(cur, query, lb, rb, errors_spent, block_id, go_right, search, blocks_length,
                                          error_left, delegate) && abort_on_hit)
        {
            return true;
        }
    }
    // Approximate search in current block.
    // i.e. blocks_length[block_id] - (rb - lb - (lb != rb)) >= min_error_left_in_block
    else if (error_left.total > 0)
    {
        // Insertion
        if (error_left.insertion > 0)
        {
            using size_type = typename cursor_t::size_type;

            size_type const lb2 = lb - !go_right;
            size_type const rb2 = rb + go_right;

            search_param error_left2{error_left};
            error_left2.total--;
            error_left2.insertion--;
            // At the end of the current block
            if (rb - lb == blocks_length[block_id])
            {
                // Leave the possibility for one or multiple deletions at the end of a block.
                // Thus do not change the direction (go_right) yet.
                // TODO: benchmark the improvement on preventing insertions followed by a deletion and vice versa. Does
                // it pay off the additional complexity and documentation for the user? (Note that the user might only
                // allow for insertions and deletion and not for mismatches).
                if (search_ss_deletion<abort_on_hit>(cur, query, lb2, rb2, errors_spent + 1, block_id, go_right, search,
                                                     blocks_length, error_left2, delegate) && abort_on_hit)
                {
                    return true;
                }
            }
            else
            {
                if (search_ss<abort_on_hit>(cur, query, lb2, rb2, errors_spent + 1, block_id, go_right, search,
                                            blocks_length, error_left2, delegate) && abort_on_hit)
                {
                    return true;
                }
            }
        }
        if (search_ss_children<abort_on_hit>(cur, query, lb, rb, errors_spent, block_id, go_right,
                                             min_error_left_in_block, search, blocks_length, error_left, delegate) &&
            abort_on_hit)
        {
            return true;
        }
    }
    return false;
}

/*!\brief Searches a query sequence in a bidirectional index using search schemes.
 * \tparam abort_on_hit     If the flag is set, the search aborts on the first hit.
 * \tparam index_t          Must model seqan3::bi_fm_index_specialisation.
 * \tparam query_t          Must model std::ranges::random_access_range over the index's alphabet.
 * \tparam search_scheme_t  Is of type `seqan3::detail::search_scheme_type` or `seqan3::detail::search_scheme_dyn_type`.
 * \tparam delegate_t       Takes `typename index_t::cursor_type` as argument.
 * \param[in] index         String index built on the text that will be searched.
 * \param[in] query         Query sequence to be searched in the index.
 * \param[in] error_left    Number of errors left for matching the remaining suffix of the query sequence.
 * \param[in] search_scheme Search scheme to be used for searching.
 * \param[in] delegate      Function that is called on every hit.
 *
 * ### Complexity
 *
 * \f$O(|query|^e)\f$ where \f$e\f$ is the total number of maximum errors.
 *
 * ### Exceptions
 *
 * Strong exception guarantee if iterating the query does not change its state and if invoking the delegate also has a
 * strong exception guarantee; basic exception guarantee otherwise.
 */
template <bool abort_on_hit, typename index_t, typename query_t, typename search_scheme_t, typename delegate_t>
inline void search_ss(index_t const & index, query_t & query, search_param const error_left,
                      search_scheme_t const & search_scheme, delegate_t && delegate)
{
    // retrieve cumulative block lengths and starting position
    auto const block_info = search_scheme_block_info(search_scheme, std::ranges::size(query));

    for (uint8_t search_id = 0; search_id < search_scheme.size(); ++search_id)
    {
        auto const & search = search_scheme[search_id];
        auto const & [blocks_length, start_pos] = block_info[search_id];

        bool const hit = search_ss<abort_on_hit>(
                             index.begin(),            // cursor on the index
                             query,                    // query to be searched
                             start_pos, start_pos + 1, // infix range already searched (open interval)
                                                       // the first character of `query` has the index 1 (not 0)
                             0,                        // errors spent
                             0,                        // current block id in search scheme
                             true,                     // search the first block from left to right
                             search, blocks_length,     // search scheme information
                             error_left,               // errors left (broken down by error types)
                             delegate                  // delegate function called on hit
                         );

        if (abort_on_hit && hit)
            return;
    }
}

/*!\brief Searches a query sequence in a bidirectional index.
 * \tparam abort_on_hit    If the flag is set, the search aborts on the first hit.
 * \tparam index_t         Must model seqan3::bi_fm_index_specialisation.
 * \tparam query_t         Must model std::ranges::random_access_range over the index's alphabet.
 * \tparam delegate_t      Takes `typename index_t::cursor_type` as argument.
 * \param[in] index        String index built on the text that will be searched.
 * \param[in] query        Query sequence to be searched in the index.
 * \param[in] error_left   Number of errors left for matching the remaining suffix of the query sequence.
 * \param[in] delegate     Function that is called on every hit.
 *
 * ### Complexity
 *
 * \f$O(|query|^e)\f$ where \f$e\f$ is the total number of maximum errors.
 *
 * ### Exceptions
 *
 * Strong exception guarantee if iterating the query does not change its state and if invoking the delegate also has a
 * strong exception guarantee; basic exception guarantee otherwise.
 */
template <bool abort_on_hit, typename index_t, typename query_t, typename delegate_t>
inline void search_algo_bi(index_t const & index, query_t & query, search_param const error_left,
                           delegate_t && delegate)
{
    switch (error_left.total)
    {
        case 0:
            search_ss<abort_on_hit>(index, query, error_left, optimum_search_scheme<0, 0>, delegate);
            break;
        case 1:
            search_ss<abort_on_hit>(index, query, error_left, optimum_search_scheme<0, 1>, delegate);
            break;
        case 2:
            search_ss<abort_on_hit>(index, query, error_left, optimum_search_scheme<0, 2>, delegate);
            break;
        case 3:
            search_ss<abort_on_hit>(index, query, error_left, optimum_search_scheme<0, 3>, delegate);
            break;
        default:
            auto const & search_scheme{compute_ss(0, error_left.total)};
            search_ss<abort_on_hit>(index, query, error_left, search_scheme, delegate);
            break;
    }
}

/*!\brief Searches a query sequence in a unidirectional index.
 *
 * \copydetails search_algo_bi
 */
template <bool abort_on_hit, typename index_t, typename query_t, typename delegate_t>
inline void search_algo_uni(index_t const & index, query_t & query, search_param const error_left,
                            delegate_t && delegate)
{
    search_trivial<abort_on_hit>(index, query, error_left, delegate);
}

/*!\brief Searches a query sequence in an index.
 *
 * \copydetails search_algo_bi
 */
template <bool abort_on_hit, typename index_t, typename query_t, typename delegate_t>
inline void search_algo(index_t const & index, query_t & query, search_param const error_left, delegate_t && delegate)
{
    if constexpr (bi_fm_index_specialisation<index_t>)
        search_algo_bi<abort_on_hit>(index, query, error_left, delegate);
    else
        search_algo_uni<abort_on_hit>(index, query, error_left, delegate);
}

//!\}

} // namespace seqan3::detail
