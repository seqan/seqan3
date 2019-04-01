// -----------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2019, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2019, Knut Reinert & MPI für molekulare Genetik
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at: https://github.com/seqan/seqan3/blob/master/LICENSE
// -----------------------------------------------------------------------------------------------------

/*!\file
 * \author Christopher Pockrandt <christopher.pockrandt AT fu-berlin.de>
 * \brief Provides the concepts for seqan3::fm_index and seqan3::bi_fm_index and its traits and cursors.
 */

#pragma once

#include <type_traits>

#include <sdsl/suffix_arrays.hpp>

#include <seqan3/core/metafunction/range.hpp>
#include <seqan3/range/concept.hpp>
#include <seqan3/range/container/concept.hpp>

namespace seqan3::detail
{

/*!\addtogroup submodule_fm_index
 * \{
 */

 // ============================================================================
 //  SdslIndex
 // ============================================================================

/*!\interface seqan3::detail::SdslIndex <>
 * \brief Concept for SDSL FM indices (which are called compressed suffix arrays in the SDSL).
 */
//!\cond
template <typename t>
SEQAN3_CONCEPT SdslIndex = requires (t sdsl_index)
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
/*!\name Requirements for seqan3::detail::SdslIndex
 * \relates seqan3::detail::SdslIndex
 * \brief The SDSL index must support the following interface to work with SeqAn3 FM indices.
 * \{
 *
 * \typedef typename t::size_type size_type
 * \brief Type for representing the size of the indexed text.
 *
 * \todo Write me.
 *
 * \}
 */

//!\}

} // namespace seqan3::detail

namespace seqan3
{

/*!\addtogroup submodule_fm_index
 * \{
 */

// ============================================================================
//  FmIndexTraits
// ============================================================================

/*!\interface seqan3::FmIndexTraits <>
 * \brief Concept for unidirectional FM Index traits.
 *
 * The traits object must contain an index type of the SDSL namespace.
 */
//!\cond
template <typename t>
SEQAN3_CONCEPT FmIndexTraits = requires (t v)
{
    typename t::sdsl_index_type;

    requires detail::SdslIndex<typename t::sdsl_index_type>;
};
//!\endcond
/*!\name Requirements for seqan3::FmIndexTraits
 * \relates seqan3::FmIndexTraits
 * \brief The SDSL index must support the following interface to work with SeqAn3.
 * \{
 *
 * \typedef typename t::sdsl_index_type sdsl_index_type
 * \brief Declares the type of the underlying SDSL index. Must satisfy the seqan3::detail::SdslIndex.
 *
 * \}
 */

// ============================================================================
//  FmIndex
// ============================================================================

/*!\interface seqan3::FmIndex <>
 * \brief Concept for unidirectional FM indices.
 *
 * This concept defines the interface for unidirectional FM indices.
 */
//!\cond
template <typename t>
SEQAN3_CONCEPT FmIndex = std::Semiregular<t> && requires (t index)
{
    typename t::text_type;
    typename t::char_type;
    typename t::size_type;
    typename t::cursor_type;

    requires sequence_container_concept<typename t::text_type>;

    // NOTE: circular dependency
    // requires FmIndexCursor<typename t::cursor_type>;

    requires requires (t index, typename t::text_type const text)
    {
        { t(text) };
        { index.construct(text) } -> void;
    };

    { index.begin() } -> typename t::cursor_type;

    { index.size()  } -> typename t::size_type;
    { index.empty() } -> bool;

    { index.load(std::string{})  } -> bool;
    { index.store(std::string{}) } -> bool;
};
//!\endcond
/*!\name Requirements for seqan3::FmIndex
 * \relates seqan3::FmIndex
 * \brief You can expect these member types and member functions on all types that satisfy seqan3::FmIndex.
 * \{
 *
 * \typedef typename t::text_type text_type
 * \brief Type of the indexed text.
 *
 * \typedef typename t::char_type char_type
 * \brief Type of the underlying character of text_type.
 *
 * \typedef typename t::size_type size_type
 * \brief Type for representing the size of the indexed text.
 *
 * \typedef typename t::cursor_type cursor_type
 * \brief Type of the unidirectional FM index cursor.
 *
 * \todo Write me!
 *
 * \}
 */

// ============================================================================
//  FmIndexCursor
// ============================================================================

/*!\interface seqan3::FmIndexCursor <>
 * \brief Concept for unidirectional FM index cursors.
 *
 * This concept defines the interface for cursors for unidirectional FM indices.
 */
//!\cond
template <typename t>
SEQAN3_CONCEPT FmIndexCursor = std::Semiregular<t> && requires (t cur)
{
    typename t::index_type;
    typename t::size_type;

    requires FmIndex<typename t::index_type>;

    requires requires (typename t::index_type const index) { { t(index) } };

    requires requires (t cur, typename t::index_type::char_type const c,
                       std::vector<typename t::index_type::char_type> const seq)
    {
        { cur.extend_right()    } -> bool;
        { cur.extend_right(c)   } -> bool;
        { cur.extend_right(seq) } -> bool;
        { cur.cycle_back()      } -> bool;
    };

    { cur.last_char()    } -> typename t::index_type::char_type;
    { cur.query_length() } -> typename t::size_type;
    { cur.query()        } -> auto;
    { *cur               } -> auto;
    { cur.count()        } -> typename t::size_type;
    { cur.locate()       } -> std::conditional_t<t::index_type::is_collection,
                                                std::vector<std::pair<typename t::size_type, typename t::size_type>>,
                                                std::vector<typename t::size_type>>;
    { cur.lazy_locate()  } -> auto;
};
//!\endcond
/*!\name Requirements for seqan3::FmIndexCursor
 * \relates seqan3::FmIndexCursor
 * \brief You can expect these member types and member functions on all types that satisfy seqan3::FmIndexCursor.
 * \{
 *
 * \typedef typename t::index_type index_type
 * \brief Type of the underlying SeqAn FM index (not the underlying SDSL index).
 *
 * \typedef typename t::size_type size_type
 * \brief Type for representing the size of the indexed text.
 *
 * \todo Write me!
 *
 * \}
 */

// ============================================================================
//  BiFmIndexTraits
// ============================================================================

/*!\interface seqan3::BiFmIndexTraits <>
 * \brief Concept for bidirectional FM Index traits.
 *
 * The traits object must contain two unidirectional FM Index traits.
 */
//!\cond
template <typename t>
SEQAN3_CONCEPT BiFmIndexTraits = requires (t v)
{
    requires FmIndexTraits<typename t::fm_index_traits>;
    requires FmIndexTraits<typename t::rev_fm_index_traits>;

    requires std::Same<typename t::fm_index_traits::sdsl_index_type::size_type,
                       typename t::rev_fm_index_traits::sdsl_index_type::size_type>;
};
//!\endcond
/*!\name Requirements for seqan3::BiFmIndexTraits
 * \relates seqan3::BiFmIndexTraits
 * \brief The bidirectional FM index traits must provide the following types:
 * \{
 *
 * \typedef typename t::fm_index_traits fm_index_traits
 * \brief Declares the type of the underlying unidirectional FM index on the original text.
 *        Must satisfy seqan3::FmIndexTraits.
 *
 * \typedef typename t::rev_fm_index_traits rev_fm_index_traits
 * \brief Declares the type of the underlying unidirectional FM index on the reversed text.
 *        Must satisfy seqan3::FmIndexTraits.
 *
 * \}
 */

// ============================================================================
//  BiFmIndex
// ============================================================================

/*!\interface seqan3::BiFmIndex <>
 * \brief Concept for bidirectional FM indices.
 *
 * This concept defines the interface for bidirectional FM indices.
 */
//!\cond
template <typename t>
SEQAN3_CONCEPT BiFmIndex = FmIndex<t> && requires (t index)
{
    typename t::cursor_type; // already required by FmIndex but has a different documentation
    typename t::fwd_cursor_type;
    typename t::rev_cursor_type;

    // NOTE: circular dependency
    // requires BiFmIndexCursor<typename t::cursor_type>;

    { index.fwd_begin() } -> typename t::fwd_cursor_type;
    { index.rev_begin() } -> typename t::rev_cursor_type;
};
//!\endcond
/*!\name Requirements for seqan3::BiFmIndex
 * \relates seqan3::BiFmIndex
 * \brief You can expect these member types and member functions on all types that satisfy seqan3::BiFmIndex.
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
 *
 * \todo Write me!
 *
 * \}
 */

// ============================================================================
//  BiFmIndexCursor
// ============================================================================

/*!\interface seqan3::BiFmIndexCursor <>
 * \brief Concept for bidirectional FM index cursors.
 *
 * This concept defines the interface for cursors for bidirectional FM indices.
 */
//!\cond
template <typename t>
SEQAN3_CONCEPT BiFmIndexCursor = FmIndexCursor<t> && requires (t cur)
{
    requires BiFmIndex<typename t::index_type>;

    requires requires (typename t::index_type const index) { { t(index) } };

    requires requires (t cur, typename t::index_type::char_type const c,
                       std::vector<typename t::index_type::char_type> const seq)
    {
        { cur.extend_left()    } -> bool;
        { cur.extend_left(c)   } -> bool;
        { cur.extend_left(seq) } -> bool;
        { cur.cycle_front()    } -> bool;
    };

};
//!\endcond
/*!\name Requirements for seqan3::BiFmIndexCursor
 * \relates seqan3::BiFmIndexCursor
 * \brief You can expect these member types and member functions on all types that satisfy
 *        seqan3::FmIndexCursor.
 * \{
 *
 * \typedef typename t::index_type index_type
 * \brief Type of the underlying SeqAn FM index (not the underlying SDSL index).
 *
 * \typedef typename t::size_type size_type
 * \brief Type for representing the size of the indexed text.
 *
 * \todo Write me!
 *
 * \}
 */

//!\}

} // namespace seqan3
