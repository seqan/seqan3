// -----------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2020, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2020, Knut Reinert & MPI für molekulare Genetik
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at: https://github.com/seqan/seqan3/blob/master/LICENSE.md
// -----------------------------------------------------------------------------------------------------

/*!\file
 * \author Christopher Pockrandt <christopher.pockrandt AT fu-berlin.de>
 * \brief Provides the internal representation of a node of the seqan3::fm_index_cursor.
 */

#pragma once

#include <tuple>
#include <type_traits>

#include <seqan3/core/platform.hpp>

namespace seqan3::detail
{

/*!\addtogroup submodule_fm_index
 * \{
 */

/*!\brief Internal representation of the node of an FM index cursor.
 * \ingroup fm_index
 * \tparam index_t The type of the underlying index; must satisfy seqan3::fm_index_specialisation.
 */
template <typename index_t>
struct fm_index_cursor_node
{
    //!\brief Type for representing positions in the indexed text.
    using size_type = typename index_t::size_type;
    /*!\brief The type of the reduced alphabet type. (The reduced alphabet might be smaller than the original alphabet
     *        in case not all possible characters occur in the indexed text.)
     */
    using sdsl_char_type = typename index_t::sdsl_char_type;

    //!\brief Left suffix array bound.
    size_type lb;
    //!\brief Right suffix array bound.
    size_type rb;
    //!\brief Depth of the node in the suffix tree, i.e. length of the searched query.
    size_type depth;
    //!\brief Label of the last edge moved down. Needed for cycle_back().
    sdsl_char_type last_char;

    //!\brief Comparison of two cursor nodes.
    bool operator==(fm_index_cursor_node const & rhs) const
    {
        // NOTE: last_char is implementation specific for cycle_back().
        // lb, rb and depth already determine the node in the suffix tree.
        // Thus there is no need to compare last_char.
        return std::tie(lb, rb, depth) == std::tie(rhs.lb, rhs.rb, rhs.depth);
    }

    //!\brief Comparison of two cursor nodes.
    bool operator!=(fm_index_cursor_node const & rhs) const
    {
        return !(*this == rhs);
    }
};

// std::tuple get_suffix_array_range(fm_index_cursor<index_t> const & it)
// {
//     return {node.lb, node.rb};
// }
//
// std::tuple get_suffix_array_range(bi_fm_index_cursor<index_t> const & it)
// {
//     return {node.lb, node.rb};
// }

//!\publicsection

//!\}

}
