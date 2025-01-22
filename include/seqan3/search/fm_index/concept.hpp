// SPDX-FileCopyrightText: 2006-2025 Knut Reinert & Freie Universität Berlin
// SPDX-FileCopyrightText: 2016-2025 Knut Reinert & MPI für molekulare Genetik
// SPDX-License-Identifier: BSD-3-Clause

/*!\file
 * \author Christopher Pockrandt <christopher.pockrandt AT fu-berlin.de>
 * \brief Provides the concept for seqan3::detail::sdsl_index.
 */

#pragma once

#include <concepts>
#include <type_traits>

#include <sdsl/suffix_arrays.hpp>

#include <seqan3/core/platform.hpp>

namespace seqan3::detail
{

// ============================================================================
//  sdsl_index
// ============================================================================

/*!\interface seqan3::detail::sdsl_index <>
 * \ingroup search_fm_index
 * \brief Concept for SDSL FM indices (which are called compressed suffix arrays in the SDSL).
 */
//!\cond
template <typename t>
concept sdsl_index = requires (t sdsl_index) {
    typename t::size_type;

    { sdsl_index.size() } -> std::same_as<typename t::size_type>;
    { sdsl_index[0] }; // suffix array access
    { sdsl_index.comp2char[0] } -> std::same_as<uint8_t>;
    { sdsl_index.char2comp[0] } -> std::same_as<uint8_t>;
    { sdsl_index.sigma };
    { sdsl_index.C[0] };

    requires requires (t sdsl_index,
                       typename t::char_type const c,
                       typename t::size_type const lb,
                       typename t::size_type const rb,
                       sdsl::int_vector<8> const text) {
        { sdsl_index.bwt.rank(lb, c) };
        { sdsl_index.wavelet_tree.lex_count(lb, rb, c) };
        { sdsl::construct_im(sdsl_index, text, 0) };
    };
};
//!\endcond
/*!\name Requirements for seqan3::detail::sdsl_index
 * \relates seqan3::detail::sdsl_index
 * \brief The SDSL index must support the following interface to work with SeqAn FM indices.
 * \{
 *
 * \typedef typename t::size_type size_type
 * \brief Type for representing the size of the indexed text.
 * \}
 */

} // namespace seqan3::detail

namespace seqan3
{
//!\brief The possible text layouts (single, collection) the seqan3::fm_index and seqan3::bi_fm_index can support.
//!\ingroup search_fm_index
enum text_layout : bool
{
    //!\brief The text is a single range.
    single,
    //!\brief The text is a range of ranges.
    collection
};

} // namespace seqan3
