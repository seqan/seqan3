// -----------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2020, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2020, Knut Reinert & MPI für molekulare Genetik
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at: https://github.com/seqan/seqan3/blob/master/LICENSE.md
// -----------------------------------------------------------------------------------------------------

/*!\file
 * \author Christopher Pockrandt <christopher.pockrandt AT fu-berlin.de>
 * \brief Provides the concepts for seqan3::fm_index and seqan3::bi_fm_index and its traits and cursors.
 */

#pragma once

#include <type_traits>

#include <sdsl/suffix_arrays.hpp>

#include <seqan3/core/type_traits/range.hpp>
#include <seqan3/range/concept.hpp>
#include <seqan3/range/container/concept.hpp>

namespace seqan3::detail
{

/*!\addtogroup submodule_fm_index
 * \{
 */

 // ============================================================================
 //  sdsl_index
 // ============================================================================

/*!\interface seqan3::detail::sdsl_index <>
 * \brief Concept for SDSL FM indices (which are called compressed suffix arrays in the SDSL).
 */
//!\cond
template <typename t>
SEQAN3_CONCEPT sdsl_index = requires (t sdsl_index)
{
    typename t::size_type;

    { sdsl_index.size() } -> typename t::size_type;
    { sdsl_index[0] }; // suffix array access
    { sdsl_index.comp2char[0] } -> uint8_t;
    { sdsl_index.char2comp[0] } -> uint8_t;
    { sdsl_index.sigma };
    { sdsl_index.C[0] };

    requires requires (t sdsl_index, typename t::char_type const c, typename t::size_type const lb,
                                     typename t::size_type const rb, sdsl::int_vector<8> const text)
    {
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

//!\}

} // namespace seqan3::detail

namespace seqan3
{

/*!\addtogroup submodule_fm_index
 * \{
 */

//!\brief The possible text layouts (single, collection) the seqan3::fm_index and seqan3::bi_fm_index can support.
enum text_layout : bool
{
    //!\brief The text is a single range.
    single,
    //!\brief The text is a range of ranges.
    collection
};

// ============================================================================
//  fm_index_specialisation
// ============================================================================

/*!\interface seqan3::fm_index_specialisation <>
 * \brief Concept for unidirectional FM indices.
 *
 * This concept defines the interface for unidirectional FM indices.
 */
//!\cond
template <typename t>
SEQAN3_CONCEPT fm_index_specialisation = std::semiregular<t> && requires (t index)
{
    typename t::alphabet_type;
    typename t::size_type;
    typename t::cursor_type;

    // NOTE: circular dependency
    // requires fm_index_cursor_specialisation<typename t::cursor_type>;
    requires requires (t index, std::conditional_t<t::text_layout_mode == text_layout::collection,
                                                   std::vector<std::vector<typename t::alphabet_type>>,
                                                   std::vector<typename t::alphabet_type>> const text)
    {
        { t(text) };
    };

    { index.cursor() } -> typename t::cursor_type;

    { index.size()  } -> typename t::size_type;
    { index.empty() } -> bool;
};
//!\endcond
/*!\name Requirements for seqan3::fm_index_specialisation
 * \relates seqan3::fm_index_specialisation
 * \brief You can expect these member types and member functions on all types that satisfy seqan3::fm_index_specialisation.
 * \{
 *
 * \typedef typename t::text_type text_type
 * \brief Type of the indexed text.
 *
 * \typedef typename t::alphabet_type alphabet_type
 * \brief Type of the underlying character of text_type.
 *
 * \typedef typename t::size_type size_type
 * \brief Type for representing the size of the indexed text.
 *
 * \typedef typename t::cursor_type cursor_type
 * \brief Type of the unidirectional FM index cursor.
 * \}
 */

// ============================================================================
//  fm_index_cursor_specialisation
// ============================================================================

/*!\interface seqan3::fm_index_cursor_specialisation <>
 * \brief Concept for unidirectional FM index cursors.
 *
 * This concept defines the interface for cursors for unidirectional FM indices.
 */
//!\cond
template <typename t>
SEQAN3_CONCEPT fm_index_cursor_specialisation = std::semiregular<t> && requires (t cur)
{
    typename t::index_type;
    typename t::size_type;

    requires fm_index_specialisation<typename t::index_type>;

    requires requires (typename t::index_type const index) { { t(index) }; };

    requires requires (t cur,
                       typename t::index_type::alphabet_type const c,
                       std::vector<typename t::index_type::alphabet_type> const seq,
                       std::conditional_t<t::index_type::text_layout_mode == text_layout::collection,
                                          std::vector<std::vector<typename t::index_type::alphabet_type>>,
                                          std::vector<typename t::index_type::alphabet_type>> const text)
    {
        { cur.extend_right()    } -> bool;
        { cur.extend_right(c)   } -> bool;
        { cur.extend_right(seq) } -> bool;
        { cur.cycle_back()      } -> bool;
        { cur.path_label(text)  } -> auto;
    };

    { cur.last_rank()    } -> typename t::size_type;
    { cur.query_length() } -> typename t::size_type;
    { cur.count()        } -> typename t::size_type;
    { cur.locate()       } -> std::conditional_t<t::index_type::text_layout_mode == text_layout::collection,
                                                std::vector<std::pair<typename t::size_type, typename t::size_type>>,
                                                std::vector<typename t::size_type>>;
    { cur.lazy_locate()  } -> auto;
};
//!\endcond
/*!\name Requirements for seqan3::fm_index_cursor_specialisation
 * \relates seqan3::fm_index_cursor_specialisation
 * \brief You can expect these member types and member functions on all types that satisfy seqan3::fm_index_cursor_specialisation.
 * \{
 *
 * \typedef typename t::index_type index_type
 * \brief Type of the underlying SeqAn FM index (not the underlying SDSL index).
 *
 * \typedef typename t::size_type size_type
 * \brief Type for representing the size of the indexed text.
 * \}
 */

// ============================================================================
//  bi_fm_index_specialisation
// ============================================================================

/*!\interface seqan3::bi_fm_index_specialisation <>
 * \brief Concept for bidirectional FM indices.
 *
 * This concept defines the interface for bidirectional FM indices.
 */
//!\cond
template <typename t>
SEQAN3_CONCEPT bi_fm_index_specialisation = fm_index_specialisation<t> && requires (t index)
{
    typename t::cursor_type; // already required by fm_index_specialisation but has a different documentation
    typename t::fwd_cursor_type;
    typename t::rev_cursor_type;

    // NOTE: circular dependency
    // requires bi_fm_index_cursor_specialisation<typename t::cursor_type>;

    { index.fwd_begin() } -> typename t::fwd_cursor_type;
    { index.rev_begin() } -> typename t::rev_cursor_type;
};
//!\endcond
/*!\name Requirements for seqan3::bi_fm_index_specialisation
 * \relates seqan3::bi_fm_index_specialisation
 * \brief You can expect these member types and member functions on all types that satisfy seqan3::bi_fm_index_specialisation.
 * \{
 *
 * \typedef typename t::cursor_type cursor_type
 * \brief Type of the bidirectional FM index cursor.
 *
 * \typedef typename t::fwd_cursor_type fwd_cursor_type
 * \brief Type of the unidirectional FM index cursor based on the unidirectional FM index on the original text.
 *
 * \typedef typename t::rev_cursor_type rev_cursor_type
 * \brief Type of the unidirectional FM index cursor based on the unidirectional FM index on the reversed text.
 * \}
 */

// ============================================================================
//  bi_fm_index_cursor_specialisation
// ============================================================================

/*!\interface seqan3::bi_fm_index_cursor_specialisation <>
 * \brief Concept for bidirectional FM index cursors.
 *
 * This concept defines the interface for cursors for bidirectional FM indices.
 */
//!\cond
template <typename t>
SEQAN3_CONCEPT bi_fm_index_cursor_specialisation = fm_index_cursor_specialisation<t> && requires (t cur)
{
    requires bi_fm_index_specialisation<typename t::index_type>;

    requires requires (typename t::index_type const index) { { t(index) }; };

    requires requires (t cur,
                       typename t::index_type::alphabet_type const c,
                       std::vector<typename t::index_type::alphabet_type> const seq)
    {
        { cur.extend_left()    } -> bool;
        { cur.extend_left(c)   } -> bool;
        { cur.extend_left(seq) } -> bool;
        { cur.cycle_front()    } -> bool;
    };

};
//!\endcond
/*!\name Requirements for seqan3::bi_fm_index_cursor_specialisation
 * \relates seqan3::bi_fm_index_cursor_specialisation
 * \brief You can expect these member types and member functions on all types that satisfy
 *        seqan3::fm_index_cursor_specialisation.
 * \{
 *
 * \typedef typename t::index_type index_type
 * \brief Type of the underlying SeqAn FM index (not the underlying SDSL index).
 *
 * \typedef typename t::size_type size_type
 * \brief Type for representing the size of the indexed text.
 * \}
 */

//!\}

} // namespace seqan3
