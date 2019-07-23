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
#include <seqan3/alignment/configuration/align_config_aligned_ends.hpp>
#include <seqan3/alignment/scoring/aminoacid_scoring_scheme.hpp>
#include <seqan3/alignment/all.hpp>
#include <seqan3/core/debug_stream.hpp>

namespace seqan3::detail
{

template <typename alphabet_t, typename TPosition>
std::vector<alphabet_t> const & getSequence(std::vector<alphabet_t> const & text, TPosition & /*pos*/)
{
    return text;
}

template <typename alphabet_t, typename T1, typename T2>
std::vector<alphabet_t> const & getSequence(std::vector<std::vector<alphabet_t>> const & text, std::pair<T1, T2> & pos)
{
//     std::cout << "Collection\n";
    return text[getSeqNo(pos)];
}

/*!\addtogroup submodule_search_algorithm
 * \{
 */

/*!\brief Searches remaining read in an text using direct comparisions between query and text or myers bitparallel algorithm.
 * \tparam abort_on_hit    If the flag is set, the search algorithm aborts on the first hit.
 * \tparam query_t         Must model std::ranges::InputRange over the index's alphabet.
 * \tparam delegate_t      Takes `index::cursor_type` as argument.
 * \param[in] cur          Cursor of atring index built on the text that will be searched.
 * \param[in] text         String that will be for the in text verification.
 * \param[in] query        Query sequence to be searched with the cursor.
 * \param[in] query_pos    Position in the query sequence indicating the prefix that has already been searched.
 * \param[in] error_left   Number of errors left for matching the remaining suffix of the query sequence.
 * \param[in] delegate     Function that is called on every hit.
 * \param[in] delegate_itv Function that is called on every hit.
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
template <bool abort_on_hit, typename text_t, typename occ_t, typename query_t, typename size_t, typename delegate_itv_t, typename align_cfg_t>
inline bool in_text_verification(text_t const & text_, occ_t & occ, query_t & query, size_t const query_pos,
                                 search_param const error_left, delegate_itv_t && delegate_itv, align_cfg_t const & align_cfg)
{
    // Select Sequence if it is a collection
    //TODO improve this
    auto const & text = getSequence(text_, occ);

    //TODO It would be nicer to check for collection but compiler does not understand that logic that pair/in of occ is always connected to collection / sequence
    //dimension_v<text_t> != 1

//         constexpr auto pos = (detail::is_type_specialisation_of<occ_t, std::pair>::value) ? occ.second : occ;
    size_t pos{};
    if constexpr (detail::is_type_specialisation_of<occ_t, std::pair>::value)
        pos = occ.second;
    else
        pos = occ;

//     std::cout << "QP: " << query_pos << "\tT:" << (int)error_left.total << "\tS: " << (int)error_left.substitution << "\tI: " << (int)error_left.insertion  << "\tD: " << (int)error_left.deletion << "\n";

    if(error_left.total <= error_left.substitution &&
        error_left.insertion == 0 && error_left.deletion == 0)
    {
//         std::cout << "Hamming\n";
        //in case comparing to text end;
        auto text_infix = text | view::slice(pos, pos + std::ranges::size(query));
        uint8_t errors = std::ranges::size(text_infix) == std::ranges::size(query) ? 0 : std::ranges::size(query) - std::ranges::size(text_infix);

        for (size_t k = query_pos; k < std::ranges::size(text_infix); ++k){
//             debug_stream << text[occ + k];
            if (text_infix[k] != query[k])
                ++errors;

            if (error_left.total < errors){
                return false;
            }
        }
        delegate_itv(occ);
        return true;
    }
    else
    {
        //Edit Distance
        //TODO remove Debug
        // Testing Alignment //TODO remove const removal as soon as slice (with gaps) works with const
//         auto text_infix = const_cast<std::vector<dna4>&>(text) | view::slice(pos + query_pos, pos + std::ranges::size(query)/* + error_left.total*/);
//         auto query_infix = const_cast<std::vector<dna4>&>(query) | view::slice(query_pos, std::ranges::size(query));

        auto text_infix = text | view::slice(pos + query_pos, pos + std::ranges::size(query) /*+ error_left.total*/);
        auto query_infix = query | view::slice(query_pos, std::ranges::size(query));

        auto results = align_pairwise(std::tie(text_infix, query_infix), align_cfg);
        auto & res = *results.begin();

            //TODO remove debug
    //         debug_stream << (text | view::slice(pos, pos + query_pos)) << "\tErrors Left:" << static_cast<int>(error_left.total) << "\n";
    //         debug_stream << (query | view::slice(0, query_pos)) << "\tError Alignment:" << (-static_cast<int>(res.score())) << "\n";
    //         debug_stream << res.alignment() << "\n";

        if(-res.score() <= static_cast<int>(error_left.total)){
    //             debug_stream << "!!!!!!!!!!!!!!!!!!!!!!!!!!Accepted!!!!!!!!!!!!!!!!!!!!!!!!!!\n";
    //             debug_stream << "ITV occ " << occ << "\n";
            delegate_itv(occ);
            return true;
        }
    return false;
    }
}

/*!\addtogroup submodule_search_algorithm
 * \{
 */

/*!\brief Searches a query sequence in an index using trivial backtracking.
 * \tparam abort_on_hit    If the flag is set, the search algorithm aborts on the first hit.
 * \tparam cursor_t        Must model seqan3::FmIndexCursor.
 * \tparam query_t         Must model std::ranges::InputRange over the index's alphabet.
 * \tparam delegate_t      Takes `index::cursor_type` as argument.
 * \param[in] cur          Cursor of atring index built on the text that will be searched.
 * \param[in] text         String that will be for the in text verification.
 * \param[in] query        Query sequence to be searched with the cursor.
 * \param[in] query_pos    Position in the query sequence indicating the prefix that has already been searched.
 * \param[in] error_left   Number of errors left for matching the remaining suffix of the query sequence.
 * \param[in] delegate     Function that is called on every hit.
 * \param[in] delegate_itv Function that is called on every hit.
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
template <bool abort_on_hit, typename query_t, typename text_t, typename cursor_t, typename delegate_t, typename delegate_itv_t, typename align_cfg_t>
inline bool search_trivial(cursor_t cur, text_t const & text, query_t & query, typename cursor_t::size_type const query_pos,
                           search_param const error_left, delegate_t && delegate, delegate_itv_t && delegate_itv, align_cfg_t const & align_cfg) noexcept(noexcept(delegate))
{
    // Exact case (end of query sequence or no errors left)
    if (query_pos == std::ranges::size(query) || (error_left.total == 0))
    {
        // If not at end of query sequence, try searching the remaining suffix without any errors.
        if (query_pos == std::ranges::size(query) || cur.extend_right(view::drop(query, query_pos)))
        {
            delegate(cur);
            return true;
        }
        return false;
    }

    // check if hamming / edit distance is applicable (cost 2-3% readlength 100 3 edit errors | for hamming distance no measureable) (gain 11% 1 Insertion 3 substitution 3 Total)(gain 10% 1I 1D 3S 3T)
    else if (query_pos >= error_left.min_step && cur.node.rb - cur.node.lb < error_left.itv_threshold && text.size() > 0 &&
        ((error_left.total <= error_left.substitution && error_left.insertion == 0 &&
        error_left.deletion == 0) || (error_left.total <= error_left.substitution &&
        error_left.total <= error_left.insertion && error_left.total <= error_left.deletion)))
    {
        for(auto occ : cur.locate())
        {
            if(in_text_verification<abort_on_hit>(text, occ, query, query_pos, error_left, delegate_itv, align_cfg))
                return true;
        }
        return false;
    }
    else
    // Approximate case
    {
        // Insertion
        if (error_left.insertion > 0)
        {
            search_param error_left2{error_left};
            error_left2.insertion--;
            error_left2.total--;

            // always perform a recursive call. Abort recursion if and only if recursive call found a hit and
            // abort_on_hit is set to true.
            if (search_trivial<abort_on_hit>(cur, text, query, query_pos + 1, error_left2, delegate, delegate_itv, align_cfg) && abort_on_hit)
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

                    if (search_trivial<abort_on_hit>(cur, text, query, query_pos + 1, error_left2, delegate, delegate_itv, align_cfg) && abort_on_hit)
                        return true;
                }

                // Deletion (Do not allow deletions at the beginning of the query sequence.)
                if (query_pos > 0)
                {
                    // Match (when error_left.substitution == 0)
                    if (error_left.substitution == 0 && cur.last_rank() == to_rank(query[query_pos]))
                    {
                        if (search_trivial<abort_on_hit>(cur, text, query, query_pos + 1, error_left, delegate, delegate_itv, align_cfg) &&
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

                        if (search_trivial<abort_on_hit>(cur, text, query, query_pos, error_left2, delegate, delegate_itv, align_cfg) && abort_on_hit)
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
                if (search_trivial<abort_on_hit>(cur, text, query, query_pos + 1, error_left, delegate, delegate_itv, align_cfg) && abort_on_hit)
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
 * \param[in] text       String that will be for the in text verification.
 * \param[in] query      Query sequence to be searched in the index.
 * \param[in] error_left Number of errors left for matching the remaining suffix of the query sequence.
 * \param[in] delegate   Function that is called on every hit.
 * \param[in] delegate_itv Function that is called on every hit.
 *
 *
 * ### Complexity
 *
 * \f$O(|query|^e)\f$ where \f$e\f$ is the maximum number of errors.
 *
 * ### Exceptions
 *
 * No-throw guarantee if invoking the delegate also guarantees no-throw.
 */
template <bool abort_on_hit, typename index_t, typename text_t, typename query_t, typename delegate_t, typename delegate_itv_t>
inline void search_trivial(index_t const & index, text_t const & text, query_t & query, search_param /*const*/ error_left,
                           delegate_t && delegate, delegate_itv_t && delegate_itv) noexcept(noexcept(delegate))
{
    // Configure the alignment kernel for in text verification.
    if constexpr (std::is_same_v<seqan3::aa27, innermost_value_type_t<text_t>> ||
                  std::is_same_v<seqan3::aa10li, innermost_value_type_t<text_t>> ||
                  std::is_same_v<seqan3::aa10murphy, innermost_value_type_t<text_t>> ||
                  std::is_same_v<seqan3::aa20, innermost_value_type_t<text_t>>){
        //myers needs to work to be competitive
        front_end_first fef{std::false_type()};
        back_end_first bef{std::true_type()};
        auto itv_cfg = align_cfg::mode{global_alignment} |
                                align_cfg::scoring{aminoacid_scoring_scheme{match_score{1}, mismatch_score{-1}}} |
                                align_cfg::gap{gap_scheme{gap_score{-1} ,gap_open_score{0}}} |
//                              align_cfg::edit | //Does not work yet
                                align_cfg::aligned_ends{end_gaps{fef, bef}}|
                                align_cfg::result{seqan3::with_score};

        //TODO turn on itv when there is myers for it hamming distance is fine.
        search_param error_left2{error_left};
        if (!(error_left.total == error_left.substitution && error_left.insertion == 0 &&
        error_left.deletion == 0) || (error_left.total <= error_left.substitution))
            error_left2.itv_threshold = 0;

        search_trivial<abort_on_hit>(index.begin(), text, query, 0, error_left2, delegate, delegate_itv, itv_cfg);

    }
    else
    {
        //TODO check edit/hamming here
        //TODO add configuration as soon as myers works with it
//         front_end_first fef{std::false_type()};
//         back_end_first bef{std::true_type()};
        auto itv_cfg  = align_cfg::edit |
    //                         align_cfg::aligned_ends{end_gaps{fef, bef}}|
                        align_cfg::result{seqan3::with_score/*with_alignment*/};

        search_trivial<abort_on_hit>(index.begin(), text, query, 0, error_left, delegate, delegate_itv, itv_cfg);
    }

}

//!\}

} // namespace seqan3::detail
