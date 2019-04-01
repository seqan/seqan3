// -----------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2019, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2019, Knut Reinert & MPI für molekulare Genetik
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at: https://github.com/seqan/seqan3/blob/master/LICENSE
// -----------------------------------------------------------------------------------------------------

/*!\file
 * \author Christopher Pockrandt <christopher.pockrandt AT fu-berlin.de>
 * \brief Provides the seqan3::bi_fm_index_cursor for searching in the bidirectional seqan3::bi_fm_index.
 */

#pragma once

#include <array>

#include <sdsl/suffix_trees.hpp>

#include <range/v3/view/slice.hpp>

#include <seqan3/alphabet/all.hpp>
#include <seqan3/core/metafunction/range.hpp>
#include <seqan3/search/fm_index/bi_fm_index.hpp>
#include <seqan3/std/ranges>

namespace seqan3
{

/*!\addtogroup submodule_fm_index
 * \{
 */

/*!\brief The SeqAn Bidirectional FM Index Cursor.
 * \implements seqan3::BiFmIndexCursor
 * \tparam index_t The type of the underlying index; must model seqan3::BiFmIndex.
 * \details
 *
 * The cursor's interface provides searching a string both from left to right as well as from right to left in the
 * indexed text. It extends the interface of the unidirectional seqan3::fm_index_cursor.
 * All methods modifying the cursor (e.g. extending by a character with extend_right()) return a `bool` value whether
 * the operation was successful or not. In case of an unsuccessful operation the cursor remains unmodified, i.e. an
 * cursor can never be in an invalid state except default constructed cursors that are always invalid.
 *
 * \if DEV
 *     The behaviour is equivalent to a prefix and suffix tree with the space and time efficiency of the underlying pure
 *     FM indices. The cursor traverses the implicit prefix and suffix trees beginning at the root node. The implicit
 *     prefix and suffix trees are not compacted, i.e. going down an edge using extend_right(char) will increase the
 *     query by only one character.
 * \endif
 *
 * The asymptotic running times for using the cursor depend on the SDSL index configuration. To determine the exact
 * running times, you have to additionally look up the running times of the used traits (configuration).
 */
template <typename index_t>
class bi_fm_index_cursor
{

public:

    //!\brief Type of the index.
    using index_type = index_t;

    /*!\name Text types
     * \{
     */
    //!\brief Type for representing positions in the indexed text.
    using size_type = typename index_type::size_type;
    //!\}

    /*!\name Cursor types
     * \{
     */
    //!\brief Type for the unidirectional cursor on the original text.
    using fwd_cursor = fm_index_cursor<fm_index<typename index_type::text_type, typename index_type::index_traits::fm_index_traits>>;
    //!\brief Type for the unidirectional cursor on the reversed text.
    using rev_cursor = fm_index_cursor<fm_index<typename index_type::rev_text_type, typename index_type::index_traits::fm_index_traits>>;
    //!\}

protected:
    //!\privatesection

    //!\brief Type of the representation of characters in the underlying SDSL index.
    using sdsl_char_type = typename index_type::sdsl_char_type;
    //!\brief Type of the alphabet size in the underlying SDSL index.
    using sdsl_sigma_type = typename index_type::sdsl_sigma_type;

    //!\brief Type of the underlying FM index.
    index_type const * index;

    /*!\name Suffix array intervals of forward and reverse cursors.
     * \{
     */
     //!\brief Left suffix array interval of the forward cursor (for extend_right).
    size_type fwd_lb;
    //!\brief Right suffix array interval of the forward cursor (for extend_right).
    size_type fwd_rb;
    //!\brief Left suffix array interval of the reverse cursor (for extend_left).
    size_type rev_lb;
    //!\brief Right suffix array interval of the reverse cursor (for extend_left).
    size_type rev_rb;
    //\}

    //!\brief Alphabet size of the index without delimiters
    sdsl_sigma_type sigma;

    //!\brief Indicates whether index is built over a collection
    static bool constexpr is_collection = dimension_v<typename index_type::text_type> == 2;

    /*!\name Information for on cycle_back() and cycle_front()
     * \brief Only stored for the cursor that has been used last to go down an edge because once one cursor is
     *        touched, the others parent information becomes invalid and cannot be used for cycle_back() anymore.
     * \{
     */

    // parent_* and _last_char only have to be stored for the (unidirectional) cursor that has been used last for
    // extend_right() or cycle_back() resp. extend_left() or cycle_front(), (i.e. either fwd or rev). Thus there is no
    // need to store it twice. Once the cursor is switched, the information becomes invalid anyway.

    //!\brief Left suffix array interval of the parent node.
    size_type parent_lb;
    //!\brief Left suffix array interval of the parent node.
    size_type parent_rb;
    //!\brief Label of the last edge moved down. Needed for cycle_back() or cycle_front().
    sdsl_char_type _last_char;
    //\}

    //!\brief Depth of the node in the suffix tree, i.e. length of the searched query.
    size_type depth; // equal for both cursors. only stored once

    // supports assertions to check whether cycle_back() resp. cycle_front() is called on the same direction as the last
    // extend_right([...]) resp. extend_left([...])
#ifndef NDEBUG
    //!\brief Stores the information which cursor has been used last for extend_*([...]) to allow for assert() in
    //        cycle_back() and cycle_front()
    bool fwd_cursor_last_used = false;
#endif

    //!\brief Helper function to recompute text positions since the indexed text is reversed.
    size_type offset() const noexcept
    {
        assert(index->size() > query_length());
        return index->size() - query_length() - 1; // since the string is reversed during construction
    }

    //!\brief Optimized bidirectional search without alphabet mapping
    template <detail::SdslIndex csa_t>
    bool bidirectional_search(csa_t const & csa, sdsl_char_type const c,
                              size_type & l_fwd, size_type & r_fwd,
                              size_type & l_bwd, size_type & r_bwd) const noexcept
    {
        assert((l_fwd <= r_fwd) && (r_fwd < csa.size()));
        assert(r_fwd + 1 >= l_fwd);
        assert(r_bwd + 1 - l_bwd == r_fwd + 1 - l_fwd);

        size_type _l_fwd, _r_fwd, _l_bwd, _r_bwd;

        size_type cc = c;
        if constexpr(!std::Same<typename csa_t::alphabet_type, sdsl::plain_byte_alphabet>)
        {
            cc = csa.char2comp[c];
            if (cc == 0 && c > 0) // [[unlikely]]
                return false;
        }

        size_type const c_begin = csa.C[cc];
        if (r_fwd + 1 - l_fwd == csa.size()) // [[unlikely]]
        {
            _l_fwd = c_begin;
            _l_bwd = c_begin;
            _r_fwd = csa.C[cc + 1] - 1;
            _r_bwd = _r_fwd;
            // if we use not the plain_byte_alphabet, we could return always return true here
        }
        else
        {
            auto      const r_s_b  = csa.wavelet_tree.lex_count(l_fwd, r_fwd + 1, c);
            size_type const rank_l = std::get<0>(r_s_b);
            size_type const s      = std::get<1>(r_s_b), b = std::get<2>(r_s_b);
            size_type const rank_r = r_fwd - l_fwd - s - b + rank_l;
            _l_fwd = c_begin + rank_l;
            _r_fwd = c_begin + rank_r;
            _l_bwd = l_bwd + s;
            _r_bwd = r_bwd - b;
        }

        if (_r_fwd >= _l_fwd)
        {
            l_fwd = _l_fwd;
            r_fwd = _r_fwd;
            l_bwd = _l_bwd;
            r_bwd = _r_bwd;
            assert(r_fwd + 1 >= l_fwd);
            assert(r_bwd + 1 - l_bwd == r_fwd + 1 - l_fwd);
            return true;
        }
        return false;
    }

    //!\brief Optimized bidirectional search for cycle_back() and cycle_front() without alphabet mapping
    template <detail::SdslIndex csa_t>
    bool bidirectional_search_cycle(csa_t const & csa, sdsl_char_type const c,
                                    size_type const l_parent, size_type const r_parent,
                                    size_type & l_fwd, size_type & r_fwd,
                                    size_type & l_bwd, size_type & r_bwd) const noexcept
    {
        assert((l_parent <= r_parent) && (r_parent < csa.size()));

        size_type c_begin;
        if constexpr(std::Same<typename csa_t::alphabet_type, sdsl::plain_byte_alphabet>)
            c_begin = csa.C[c]; // TODO: check whether this can be removed
        else
            c_begin = csa.C[csa.char2comp[c]];

        auto const r_s_b = csa.wavelet_tree.lex_count(l_parent, r_parent + 1, c);
        size_type const s      = std::get<1>(r_s_b),
                        b      = std::get<2>(r_s_b),
                        rank_l = std::get<0>(r_s_b),
                        rank_r = r_parent - l_parent - s - b + rank_l;

        size_type const _l_fwd = c_begin + rank_l;
        size_type const _r_fwd = c_begin + rank_r;
        size_type const _l_bwd = r_bwd + 1;
        size_type const _r_bwd = r_bwd + 1 + rank_r - rank_l;

        if (_r_fwd >= _l_fwd)
        {
            l_fwd = _l_fwd;
            r_fwd = _r_fwd;
            l_bwd = _l_bwd;
            r_bwd = _r_bwd;
            assert(r_fwd + 1 >= l_fwd);
            assert(r_bwd + 1 - l_bwd == r_fwd + 1 - l_fwd);
            return true;
        }
        return false;
    }

public:

    /*!\name Constructors, destructor and assignment
     * \{
     */
    //!\brief Default constructor. Accessing member functions on a default constructed object is undefined behavior.
    //        Default construction is necessary to make this class semi-regular and e.g., to allow construction of
    //        std::array of cursors.
    bi_fm_index_cursor() noexcept = default;                                       //!< Default constructor.
    bi_fm_index_cursor(bi_fm_index_cursor const &) noexcept = default;             //!< Copy constructor.
    bi_fm_index_cursor & operator=(bi_fm_index_cursor const &) noexcept = default; //!< Copy assignment.
    bi_fm_index_cursor(bi_fm_index_cursor &&) noexcept = default;                  //!< Move constructor.
    bi_fm_index_cursor & operator=(bi_fm_index_cursor &&) noexcept = default;      //!< Move assignment.
    ~bi_fm_index_cursor() = default;                                               //!< Destructor.

    //! \brief Construct from given index.
    bi_fm_index_cursor(index_t const & _index) noexcept : index(&_index),
                                                          fwd_lb(0), fwd_rb(_index.size() - 1),
                                                          rev_lb(0), rev_rb(_index.size() - 1),
                                                          sigma(_index.fwd_fm.index.sigma - is_collection),
                                                          depth(0)
    {}
    //\}

    /*!\brief Compares two cursors.
     * \param[in] rhs Other cursor to compare it to.
     * \returns `true` if both cursors are equal, `false` otherwise.
     *
     * ### Complexity
     *
     * Constant.
     *
     * ### Exceptions
     *
     * No-throw guarantee.
     */
    bool operator==(bi_fm_index_cursor const & rhs) const noexcept
    {
        assert(index != nullptr);
        // equal SA interval implies equal parent node information (or both are root nodes)
        assert(!(fwd_lb == rhs.fwd_lb && fwd_rb == rhs.fwd_rb && depth == rhs.depth) ||
               (depth == 0) ||
               (parent_lb == rhs.parent_lb && parent_rb == rhs.parent_rb && _last_char == rhs._last_char));

        return std::tie(fwd_lb, fwd_rb, depth) == std::tie(rhs.fwd_lb, rhs.fwd_rb, rhs.depth);
    }

    /*!\brief Compares two cursors.
     * \param[in] rhs Other cursor to compare it to.
     * \returns `true` if the cursors are not equal, `false` otherwise.
     *
     * ### Complexity
     *
     * Constant.
     *
     * ### Exceptions
     *
     * No-throw guarantee.
     */
    bool operator!=(bi_fm_index_cursor const & rhs) const noexcept
    {
        assert(index != nullptr);

        return !(*this == rhs);
    }

    /*!\brief Tries to extend the query by the smallest possible character to the right such that the query is found in
     *        the text.
     *        \if DEV
     *            Goes down the leftmost (i.e. lexicographically smallest) edge.
     *        \endif
     * \returns `true` if the cursor could extend the query successfully.
     *
     * ### Complexity
     *
     * \f$O(\Sigma) * O(T_{BACKWARD\_SEARCH})\f$
     *
     * It scans linearly over the alphabet until it finds the smallest character that is represented by an edge.
     *
     * ### Exceptions
     *
     * No-throw guarantee.
     */
    bool extend_right() noexcept
    {
    #ifndef NDEBUG
        fwd_cursor_last_used = true;
    #endif

        assert(index != nullptr);

        size_type new_parent_lb = fwd_lb, new_parent_rb = fwd_rb;

        sdsl_char_type c = 1; // NOTE: start with 0 or 1 depending on implicit_sentintel
        while (c < sigma &&
               !bidirectional_search(index->fwd_fm.index, index->fwd_fm.index.comp2char[c],
                                     fwd_lb, fwd_rb, rev_lb, rev_rb))
        {
            ++c;
        }

        if (c != sigma)
        {
            parent_lb = new_parent_lb;
            parent_rb = new_parent_rb;

            _last_char = c;
            ++depth;

            return true;
        }
        return false;
    }

    /*!\brief Tries to extend the query by the smallest possible character to the left such that the query is found in
     *        the text.
     *        \if DEV
     *            Goes down the leftmost (i.e. lexicographically smallest) edge in the reverse cursor.
     *        \endif
     * \returns `true` if the cursor could extend the query successfully.
     *
     * ### Complexity
     *
     * \f$O(\Sigma) * O(T_{BACKWARD\_SEARCH})\f$
     *
     * It scans linearly over the alphabet until it finds the smallest character that is represented by an edge.
     *
     * ### Exceptions
     *
     * No-throw guarantee.
     */
    bool extend_left() noexcept
    {
    #ifndef NDEBUG
        fwd_cursor_last_used = false;
    #endif

        assert(index != nullptr);

        size_type new_parent_lb = rev_lb, new_parent_rb = rev_rb;

        sdsl_char_type c = 1; // NOTE: start with 0 or 1 depending on implicit_sentintel
        while (c < sigma &&
               !bidirectional_search(index->rev_fm.index, index->rev_fm.index.comp2char[c],
                                     rev_lb, rev_rb, fwd_lb, fwd_rb))
        {
            ++c;
        }

        if (c != sigma)
        {
            parent_lb = new_parent_lb;
            parent_rb = new_parent_rb;

            _last_char = c;
            ++depth;

            return true;
        }
        return false;
    }

    /*!\brief Tries to extend the query by the character `c` to the right.
     * \tparam char_t Type of the character needs to be convertible to the character type `char_type` of the indexed
     *                text.
     * \param[in] c Character to extend the query with to the right.
     * \returns `true` if the cursor could extend the query successfully.
     *
     * ### Complexity
     *
     * \f$O(T_{BACKWARD\_SEARCH})\f$
     *
     * ### Exceptions
     *
     * No-throw guarantee.
     */
    template <Alphabet char_t>
    //!\cond
        requires ImplicitlyConvertibleTo<char_t, typename index_t::char_type>
    //!\endcond
    bool extend_right(char_t const c) noexcept
    {
    #ifndef NDEBUG
        fwd_cursor_last_used = true;
    #endif

        assert(index != nullptr);

        size_type new_parent_lb = fwd_lb, new_parent_rb = fwd_rb;

        auto c_char = to_rank(c) + 1;
        if (bidirectional_search(index->fwd_fm.index, c_char, fwd_lb, fwd_rb, rev_lb, rev_rb))
        {
            parent_lb = new_parent_lb;
            parent_rb = new_parent_rb;

            _last_char = c_char;
            ++depth;

            return true;
        }
        return false;
    }

    /*!\brief Tries to extend the query by the character `c` to the left.
     * \tparam char_t Type of the character needs to be convertible to the character type `char_type` of the indexed
     *                text.
     * \param[in] c Character to extend the query with to the left.
     * \returns `true` if the cursor could extend the query successfully.
     *
     * ### Complexity
     *
     * \f$O(T_{BACKWARD\_SEARCH})\f$
     *
     * ### Exceptions
     *
     * No-throw guarantee.
     */
    template <Alphabet char_t>
    //!\cond
       requires ImplicitlyConvertibleTo<char_t, typename index_t::char_type>
    //!\endcond
    bool extend_left(char_t const c) noexcept
    {
    #ifndef NDEBUG
        fwd_cursor_last_used = false;
    #endif

        assert(index != nullptr);

        size_type new_parent_lb = rev_lb, new_parent_rb = rev_rb;

        auto c_char = to_rank(c) + 1;
        if (bidirectional_search(index->rev_fm.index, c_char, rev_lb, rev_rb, fwd_lb, fwd_rb))
        {
            parent_lb = new_parent_lb;
            parent_rb = new_parent_rb;

            _last_char = c_char;
            ++depth;

            return true;
        }
        return false;
    }

    /*!\brief Tries to extend the query by `seq` to the right.
     * \tparam seq_t The type of range of the sequence to search; must model std::ranges::RandomAccessRange.
     * \param[in] seq Sequence to extend the query with to the right.
     * \returns `true` if the cursor could extend the query successfully.
     *
     * If extending fails in the middle of the sequence, all previous computations are rewound to restore the cursor's
     * state before calling this method.
     *
     * ### Complexity
     *
     * \f$|seq| * O(T_{BACKWARD\_SEARCH})\f$
     *
     * ### Exceptions
     *
     * No-throw guarantee.
     */
    template <std::ranges::RandomAccessRange seq_t>
    //!\cond
        requires ImplicitlyConvertibleTo<innermost_value_type_t<seq_t>, typename index_t::char_type>
    //!\endcond
    bool extend_right(seq_t && seq) noexcept
    {
        assert(index != nullptr);

        auto first = std::ranges::begin(seq);
        auto last = std::ranges::end(seq);

    #ifndef NDEBUG
        fwd_cursor_last_used = (first != last); // only if seq was not empty
    #endif

        size_type _fwd_lb = fwd_lb, _fwd_rb = fwd_rb, _rev_lb = rev_lb, _rev_rb = rev_rb;
        size_type new_parent_lb = parent_lb, new_parent_rb = parent_rb;
        sdsl_char_type c = _last_char;

        for (auto it = first; it != last; ++it)
        {
            c = to_rank(*it) + 1;

            new_parent_lb = _fwd_lb;
            new_parent_rb = _fwd_rb;
            if (!bidirectional_search(index->fwd_fm.index, c, _fwd_lb, _fwd_rb, _rev_lb, _rev_rb))
                return false;
        }

        fwd_lb = _fwd_lb;
        fwd_rb = _fwd_rb;
        rev_lb = _rev_lb;
        rev_rb = _rev_rb;

        parent_lb = new_parent_lb;
        parent_rb = new_parent_rb;

        _last_char = c;
        depth += last - first;

        return true;
    }

    /*!\brief Tries to extend the query by `seq` to the left.
     * \tparam seq_t The type of range of the sequence to search; must model std::ranges::RandomAccessRange.
     * \param[in] seq Sequence to extend the query with to the left (starting from right to left, see example).
     * \returns `true` if the cursor could extend the query successfully.
     *
     * If extending fails in the middle of the sequence, all previous computations are rewound to restore the cursor's
     * state before calling this method.
     *
     * Example:
     *
     * \snippet test/snippet/search/bi_fm_index_cursor.cpp extend_left_seq
     *
     * ### Complexity
     *
     * \f$|seq| * O(T_{BACKWARD\_SEARCH})\f$
     *
     * ### Exceptions
     *
     * No-throw guarantee.
     */
    template <std::ranges::RandomAccessRange seq_t>
    //!\cond
       requires ImplicitlyConvertibleTo<innermost_value_type_t<seq_t>, typename index_t::char_type>
    //!\endcond
    bool extend_left(seq_t && seq) noexcept
    {
        assert(index != nullptr);

        auto rev_seq = std::view::reverse(seq);
        auto first = std::ranges::begin(rev_seq);
        auto last = std::ranges::end(rev_seq);

    #ifndef NDEBUG
        if (first != last) // only if seq was not empty
            fwd_cursor_last_used = false;
    #endif

        size_type _fwd_lb = fwd_lb, _fwd_rb = fwd_rb,
                  _rev_lb = rev_lb, _rev_rb = rev_rb;
        size_type new_parent_lb = parent_lb, new_parent_rb = parent_rb;
        sdsl_char_type c = _last_char;

        for (auto it = first; it != last; ++it)
        {
            c = to_rank(*it) + 1;

            new_parent_lb = _rev_lb;
            new_parent_rb = _rev_rb;
            if (!bidirectional_search(index->rev_fm.index, c, _rev_lb, _rev_rb, _fwd_lb, _fwd_rb))
                return false;
        }

        fwd_lb = _fwd_lb;
        fwd_rb = _fwd_rb;
        rev_lb = _rev_lb;
        rev_rb = _rev_rb;

        parent_lb = new_parent_lb;
        parent_rb = new_parent_rb;
        _last_char = c;
        depth += last - first;

        return true;
    }

    /*!\brief Tries to replace the rightmost character of the query by the next lexicographically larger character such
     *        that the query is found in the text.
     *        \if DEV
     *            Moves the cursor to the right sibling of the current suffix tree node. It would be equivalent to
     *            going up an edge and going down that edge with the smallest character that is larger than the
     *            previous searched character. Calling cycle_*() on an cursor pointing to the root node is undefined
     *            behaviour!
     *        \endif
     * \returns `true` if there exists a query in the text where the rightmost character of the query is
     *          lexicographically larger than the current rightmost character of the query.
     *
     * Example:
     *
     * \snippet test/snippet/search/bi_fm_index_cursor.cpp cycle
     *
     * ### Complexity
     *
     * \f$O(\Sigma) * O(T_{BACKWARD\_SEARCH})\f$
     *
     * It scans linearly over the alphabet starting from the rightmost character until it finds the query with a larger
     * rightmost character.
     *
     * ### Exceptions
     *
     * No-throw guarantee.
     */
    bool cycle_back() noexcept
    {
    #ifndef NDEBUG
        // cycle_back() can only be used if the last extension was to the right.
        assert(fwd_cursor_last_used);
    #endif

        assert(index != nullptr && query_length() > 0);

        sdsl_char_type c = _last_char + 1;

        while (c < sigma &&
               !bidirectional_search_cycle(index->fwd_fm.index, index->fwd_fm.index.comp2char[c],
                                           parent_lb, parent_rb, fwd_lb, fwd_rb, rev_lb, rev_rb))
        {
            ++c;
        }

        if (c != sigma)
        {
            _last_char = c;

            return true;
        }
        return false;
    }

    /*!\brief Tries to replace the leftmost character of the query by the next lexicographically larger character such
     *        that the query is found in the text.
     *        \if DEV
     *            Moves the cursor to the right sibling of the current suffix tree node. It would be equivalent to
     *            going up an edge and going down that edge with the smallest character that is larger than the
     *            previous searched character. Calling cycle_*() on an cursor pointing to the root node is undefined
     *            behaviour!
     *        \endif
     * \returns `true` if there exists a query in the text where the leftmost character of the query is
     *          lexicographically larger than the current leftmost character of the query.
     *
     * Example:
     *
     * \snippet test/snippet/search/bi_fm_index_cursor.cpp cycle
     *
     * ### Complexity
     *
     * \f$O(\Sigma) * O(T_{BACKWARD\_SEARCH})\f$
     *
     * It scans linearly over the alphabet starting from the leftmost character until it finds the query with a larger
     * leftmost character.
     *
     * ### Exceptions
     *
     * No-throw guarantee.
     */
    bool cycle_front() noexcept
    {
    #ifndef NDEBUG
        // cycle_front() can only be used if the last extension was to the left.
        assert(!fwd_cursor_last_used);
    #endif

        assert(index != nullptr && query_length() > 0);

        sdsl_char_type c = _last_char + 1;
        while (c < sigma &&
               !bidirectional_search_cycle(index->rev_fm.index, index->rev_fm.index.comp2char[c],
                                           parent_lb, parent_rb, rev_lb, rev_rb, fwd_lb, fwd_rb))
        {
            ++c;
        }

        if (c != sigma)
        {
            _last_char = c;

            return true;
        }
        return false;
    }


    /*!\brief Outputs the rightmost respectively leftmost character depending on whether extend_right() or extend_left()
     *        has been called last.
     * \returns Rightmost or leftmost character.
     *
     * Example:
     *
     * \snippet test/snippet/search/bi_fm_index_cursor.cpp cycle
     *
     * ### Complexity
     *
     * Constant.
     *
     * ### Exceptions
     *
     * No-throw guarantee.
     */
    typename index_t::char_type last_char() noexcept
    {
        assert(index != nullptr && query_length() > 0);

        typename index_t::char_type c;
        assign_rank(c, index->fwd_fm.index.comp2char[_last_char] - 1); // text is not allowed to contain ranks of 0
        return c;
    }

    /*!\brief Returns the depth of the cursor node in the implicit suffix tree, i.e. the length of the sequence
     *        searched.
     * \returns Length of searched sequence.
     *
     * ### Complexity
     *
     * Constant.
     *
     * ### Exceptions
     *
     * No-throw guarantee.
     */
    size_type query_length() const noexcept
    {
        assert(index != nullptr);
        // depth == 0 -> root node
        assert(depth != 0 ||
               (fwd_lb == rev_lb &&
                fwd_rb == rev_rb &&
                fwd_lb == 0 &&
                fwd_rb == index->size() - 1));

        return depth;
    }

    /*!\brief Returns a unidirectional seqan3::fm_index_cursor on the original text. query() on the returned
     *        unidirectional index cursor will be equal to query() on the bidirectional index cursor.
     *        cycle_back() and last_char() will be undefined behavior if the last extension on the bidirectional
     *        FM index has been to the left. The behavior will be well-defined after the first extension to the right
     *        on the unidirectional index.
     * \returns Returns a unidirectional seqan3::fm_index_cursor on the index of the original text.
     *
     * Example:
     *
     * \snippet test/snippet/search/bi_fm_index_cursor.cpp to_fwd_cursor
     *
     * ### Complexity
     *
     * Constant.
     *
     * ### Exceptions
     *
     * No-throw guarantee.
     */
    fwd_cursor to_fwd_cursor() const noexcept
    {
        assert(index != nullptr);

        fwd_cursor cur{index->fwd_fm};
        cur.parent_lb = parent_lb;
        cur.parent_rb = parent_rb;
        cur.node = {fwd_lb, fwd_rb, depth, _last_char};

    #ifndef NDEBUG
        if (!fwd_cursor_last_used)
        {
            // invalidate parent interval
            cur.parent_lb = 1;
            cur.parent_rb = 0;
        }
    #endif

        return cur;
    }

    /*!\brief Returns a unidirectional seqan3::fm_index_cursor on the reversed text. query() on the returned
     *        unidirectional index cursor will be equal to reversing query() on the bidirectional index cursor.
     *        Note that because of the text being reversed, extend_right() resp. cycle_back()
     *        correspond to extend_left() resp. cycle_front() on the bidirectional index cursor.
     *        Furthermore cycle_back() and last_char() will be undefined behavior if the last extension on the
     *        bidirectional FM index has been to the left. The behavior will be well-defined after the first
     *        extension to the right on the unidirectional index.
     * \returns Returns a unidirectional seqan3::fm_index_cursor on the index of the reversed text.
     *
     * Example:
     *
     * \snippet test/snippet/search/bi_fm_index_cursor.cpp to_rev_cursor
     *
     * \attention When the index is built for text collections, the returned text IDs will be reversed.
     *
     * \snippet test/snippet/search/bi_fm_index_cursor.cpp to_rev_cursor_collection
     * ### Complexity
     *
     * Constant.
     *
     * ### Exceptions
     *
     * No-throw guarantee.
     */
    rev_cursor to_rev_cursor() const noexcept
    {
        assert(index != nullptr);

        rev_cursor cur{index->rev_fm};
        cur.parent_lb = parent_lb;
        cur.parent_rb = parent_rb;
        cur.node = {rev_lb, rev_rb, depth, _last_char};

    #ifndef NDEBUG
        if (fwd_cursor_last_used)
        {
            // invalidate parent interval
            cur.parent_lb = 1;
            cur.parent_rb = 0;
        }
    #endif

        return cur;
    }

    /*!\brief Returns the searched query.
     *        \if DEV
     *            Returns the concatenation of all edges from the root node to the cursors current node.
     *        \endif
     *
     * ### Complexity
     *
     * \f$O(SAMPLING\_RATE * T_{BACKWARD\_SEARCH}) + query\_length()\f$
     *
     * ### Exceptions
     *
     * No-throw guarantee.
     */
    auto query() const noexcept
    //!\cond
        requires !is_collection
    //!\endcond
    {
        assert(index != nullptr && index->text != nullptr);

        size_type const query_begin = offset() - index->fwd_fm.index[fwd_lb];
        return *index->text | ranges::view::slice(query_begin, query_begin + query_length());
    }

    //!\overload
    auto query() const noexcept
    //!\cond
        requires is_collection
    //!\endcond
    {
        assert(index != nullptr && index->text != nullptr);

        size_type const loc = offset() - index->fwd_fm.index[fwd_lb];
        size_type const query_begin = loc - index->fwd_fm.text_begin_rs.rank(loc + 1) + 1; // Substract delimiters
        return *index->text | std::view::join | ranges::view::slice(query_begin, query_begin + query_length());
    }

    //!\copydoc query()
    auto operator*() const noexcept
    {
        assert(index != nullptr && index->text != nullptr);

        return query();
    }

    /*!\brief Counts the number of occurrences of the searched query in the text.
     * \returns Number of occurrences of the searched query in the text.
     *
     * ### Complexity
     *
     * Constant.
     *
     * ### Exceptions
     *
     * No-throw guarantee.
     */
    size_type count() const noexcept
    {
        assert(index != nullptr && (1 + fwd_rb - fwd_lb == 1 + rev_rb - rev_lb));

        return 1 + fwd_rb - fwd_lb;
    }

    /*!\brief Locates the occurrences of the searched query in the text.
     * \returns Positions in the text.
     *
     * ### Complexity
     *
     * \f$count() * O(T_{BACKWARD\_SEARCH} * SAMPLING\_RATE)\f$
     *
     * ### Exceptions
     *
     * Strong exception guarantee (no data is modified in case an exception is thrown).
     */
    std::vector<size_type> locate() const
    //!\cond
        requires !is_collection
    //!\endcond
    {
        assert(index != nullptr);

        std::vector<size_type> occ(count());
        for (size_type i = 0; i < occ.size(); ++i)
        {
            occ[i] = offset() - index->fwd_fm.index[fwd_lb + i];
        }
        return occ;
    }

    //!\overload
    std::vector<std::pair<size_type, size_type>> locate() const
    //!\cond
        requires is_collection
    //!\endcond
    {
        assert(index != nullptr);

        std::vector<std::pair<size_type, size_type>> occ;
        occ.reserve(count());
        for (size_type i = 0; i < count(); ++i)
        {
            size_type loc = offset() - index->fwd_fm.index[fwd_lb + i];
            size_type sequence_rank = index->fwd_fm.text_begin_rs.rank(loc + 1);
            size_type sequence_position = loc - index->fwd_fm.text_begin_ss.select(sequence_rank);
            occ.emplace_back(sequence_rank - 1, sequence_position);
        }
        return occ;
    }

    /*!\brief Locates the occurrences of the searched query in the text on demand, i.e. a ranges::view is returned
     *        and every position is located once it is accessed.
     * \returns Positions in the text.
     *
     * ### Complexity
     *
     * \f$count() * O(T_{BACKWARD\_SEARCH} * SAMPLING\_RATE)\f$
     *
     * ### Exceptions
     *
     * Strong exception guarantee (no data is modified in case an exception is thrown).
     */
    auto lazy_locate() const
    //!\cond
        requires !is_collection
    //!\endcond
    {
        assert(index != nullptr);

        return std::view::iota(fwd_lb, fwd_lb + count())
             | std::view::transform([*this, _offset = offset()] (auto sa_pos)
               {
                   return _offset - index->fwd_fm.index[sa_pos];
               });
    }

    //!\overload
    auto lazy_locate() const
    //!\cond
        requires is_collection
    //!\endcond
    {
        assert(index != nullptr);

        return std::view::iota(fwd_lb, fwd_lb + count())
               | std::view::transform([*this, _offset = offset()] (auto sa_pos)
               {
                   return _offset - index->fwd_fm.index[sa_pos];
               })
               | std::view::transform([*this] (auto loc)
               {
                   size_type sequence_rank = index->fwd_fm.text_begin_rs.rank(loc + 1);
                   size_type sequence_position = loc - index->fwd_fm.text_begin_ss.select(sequence_rank);
                   return std::make_pair(sequence_rank-1, sequence_position);
               });
    }

};

//!\}

} // namespace seqan3
