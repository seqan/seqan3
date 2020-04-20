// -----------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2020, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2020, Knut Reinert & MPI für molekulare Genetik
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at: https://github.com/seqan/seqan3/blob/master/LICENSE.md
// -----------------------------------------------------------------------------------------------------

/*!\file
 * \author Christopher Pockrandt <christopher.pockrandt AT fu-berlin.de>
 * \brief Provides the seqan3::fm_index_cursor for searching in the unidirectional seqan3::fm_index.
 */

#pragma once

#include <array>
#include <type_traits>

#include <sdsl/suffix_trees.hpp>

#include <range/v3/view/slice.hpp>

#include <seqan3/alphabet/adaptation/char.hpp>
#include <seqan3/alphabet/concept.hpp>
#include <seqan3/core/type_traits/range.hpp>
#include <seqan3/range/views/slice.hpp>
#include <seqan3/search/fm_index/detail/fm_index_cursor.hpp>
#include <seqan3/search/fm_index/fm_index.hpp>
#include <seqan3/std/ranges>

namespace seqan3
{

/*!\addtogroup submodule_fm_index
 * \{
 */

/*!\brief The SeqAn FM Index Cursor.
 * \implements seqan3::fm_index_cursor_specialisation
 * \tparam index_t The type of the underlying index; must model seqan3::fm_index_specialisation.
 * \details
 *
 * The cursor's interface provides searching a string from left to right in the indexed text.
 * All methods modifying the cursor (e.g. extending by a character with extend_right()) return a `bool` value whether
 * the operation was successful or not. In case of an unsuccessful operation the cursor remains unmodified, i.e. an
 * cursor can never be in an invalid state except default constructed cursors that are always invalid.
 *
 * \if DEV
 *     The behaviour is equivalent to a suffix tree with the space and time efficiency of the underlying pure FM index.
 *     The cursor traverses the implicit suffix tree beginning at the root node.
 *     The implicit suffix tree is not compacted, i.e. going down an edge using extend_right(char) will increase the
 *     query by only one character.
 * \endif
 *
 * The asymptotic running times for using the cursor depend on the SDSL index configuration. To determine the exact
 * running times, you have to additionally look up the running times of the used traits (configuration).
 */
template <typename index_t>
class fm_index_cursor
{
public:

    /*!\name Member types
     * \{
     */
    //!\brief Type of the index.
    using index_type = index_t;
    //!\brief Type for representing positions in the indexed text.
    using size_type = typename index_type::size_type;
    //!\}

private:
    //!\privatesection

    /*!\name Member types
     * \{
     */
    //!\brief Type of the representation of a suffix tree node.
    using node_type = detail::fm_index_cursor_node<index_t>;
    //!\brief Type of the representation of characters in the underlying SDSL index.
    using sdsl_char_type = typename index_type::sdsl_char_type;
    //!\brief Type of the SDSL index.
    using sdsl_index_type = typename index_t::sdsl_index_type;
    //!\brief Type of the alphabet size in the underlying SDSL index.
    using sdsl_sigma_type = typename index_type::sdsl_sigma_type;
    //!\brief Alphabet type of the index.
    using index_alphabet_type = typename index_t::alphabet_type;
    //!\}

    //!\brief Underlying FM index.
    index_type const * index{nullptr};
    //!\brief Left suffix array interval of the parent node. Needed for cycle_back().
    size_type parent_lb{};
    //!\brief Right suffix array interval of the parent node. Needed for cycle_back().
    size_type parent_rb{};
    //!\brief Underlying index from the SDSL.
    node_type node{};
    //!\brief Alphabet size of the index without delimiters
    sdsl_sigma_type sigma{};

    template <typename _index_t>
    friend class bi_fm_index_cursor;

    //!\brief Helper function to recompute text positions since the indexed text is reversed.
    size_type offset() const noexcept
    {
        assert(index->index.size() > query_length());
        return index->index.size() - query_length() - 1; // since the string is reversed during construction
    }

    //!\brief Optimized backward search without alphabet mapping
    bool backward_search(sdsl_index_type const & csa,
                         sdsl_char_type const c,
                         size_type & l,
                         size_type & r) const noexcept
    {
        assert(l <= r && r < csa.size());

        size_type _l, _r;

        size_type cc = c;
        if constexpr(!std::same_as<index_alphabet_type, sdsl::plain_byte_alphabet>)
        {
            cc = csa.char2comp[c];
            if (cc == 0 && c > 0) // [[unlikely]]
                return false;
        }

        size_type const c_begin = csa.C[cc];
        if (l == 0 && r + 1 == csa.size()) // [[unlikely]]
        {
            _l = c_begin;
            _r = csa.C[cc + 1] - 1;
            // if we use not the plain_byte_alphabet, we could return always return true here
        }
        else
        {
            _l = c_begin + csa.bwt.rank(l, c);		   // count c in bwt[0..l-1]
            _r = c_begin + csa.bwt.rank(r + 1, c) - 1; // count c in bwt[0..r]
        }

        if (_r >= _l)
        {
            r = _r;
            l = _l;
            assert(r + 1 - l >= 0);
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
    //        std::array of iterators.
    fm_index_cursor() noexcept = default;                                    //!< Defaulted.
    fm_index_cursor(fm_index_cursor const &) noexcept = default;             //!< Defaulted.
    fm_index_cursor & operator=(fm_index_cursor const &) noexcept = default; //!< Defaulted.
    fm_index_cursor(fm_index_cursor &&) noexcept = default;                  //!< Defaulted.
    fm_index_cursor & operator=(fm_index_cursor &&) noexcept = default;      //!< Defaulted.
    ~fm_index_cursor() = default;                                            //!< Defaulted.

    //! \brief Construct from given index.
    fm_index_cursor(index_t const & _index) noexcept :
        index(&_index),
        node({0, _index.index.size() - 1, 0, 0}),
        sigma(_index.index.sigma - index_t::text_layout_mode)
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
    bool operator==(fm_index_cursor const & rhs) const noexcept
    {
        assert(index != nullptr);
        assert(node != rhs.node || (query_length() == 0 || (parent_lb == rhs.parent_lb && parent_rb == rhs.parent_rb)));

        // position in the implicit suffix tree is defined by the SA interval and depth.
        // No need to compare parent intervals
        return node == rhs.node;
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
    bool operator!=(fm_index_cursor const & rhs) const noexcept
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
        // TODO: specialize extend_right() and cycle_back() for EPR-dictionaries
        // store all cursors at once in a private std::array of cursors
        assert(index != nullptr);

        sdsl_char_type c = 1; // NOTE: start with 0 or 1 depending on implicit_sentintel
        size_type _lb = node.lb, _rb = node.rb;
        while (c < sigma && !backward_search(index->index, index->index.comp2char[c], _lb, _rb))
        {
            ++c;
        }

        if (c != sigma)
        {
            parent_lb = node.lb;
            parent_rb = node.rb;
            node = {_lb, _rb, node.depth + 1, c};
            return true;
        }
        return false;
    }

    /*!\brief Tries to extend the query by the character `c` to the right.
     * \tparam char_t Type of the character needs to be convertible to the character type `char_type` of the index.
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
    template <typename char_t>
    bool extend_right(char_t const c) noexcept
    {
        static_assert(std::convertible_to<char_t, index_alphabet_type>,
                     "The character must be convertible to the alphabet of the index.");

        assert(index != nullptr);
        // The rank cannot exceed 255 for single text and 254 for text collections as they are reserved as sentinels
        // for the indexed text.
        assert(seqan3::to_rank(static_cast<index_alphabet_type>(c)) <
               ((index_type::text_layout_mode == text_layout::single) ? 255 : 254));

        size_type _lb = node.lb, _rb = node.rb;

        sdsl_char_type c_char = seqan3::to_rank(static_cast<index_alphabet_type>(c)) + 1;

        if (backward_search(index->index, c_char, _lb, _rb))
        {
            parent_lb = node.lb;
            parent_rb = node.rb;
            node = {_lb, _rb, node.depth + 1, c_char};
            return true;
        }
        return false;
    }

    //!\overload
    template <typename char_type>
    //!\cond
        requires detail::is_char_adaptation_v<char_type>
    //!\endcond
    bool extend_right(char_type const * cstring) noexcept
    {
        return extend_right(std::basic_string_view<char_type>{cstring});
    }

    /*!\brief Tries to extend the query by `seq` to the right.
     * \tparam seq_t The type of range of the sequence to search; must model std::ranges::forward_range.
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
    template <std::ranges::range seq_t>
    bool extend_right(seq_t && seq) noexcept
    {
        static_assert(std::ranges::forward_range<seq_t>, "The query must model forward_range.");
        static_assert(std::convertible_to<innermost_value_type_t<seq_t>, index_alphabet_type>,
                     "The alphabet of the sequence must be convertible to the alphabet of the index.");

        assert(index != nullptr); // range must not be empty!

        size_type _lb = node.lb, _rb = node.rb;
        size_type new_parent_lb = parent_lb, new_parent_rb = parent_rb;

        sdsl_char_type c{};
        size_t len{0};

        for (auto it = std::ranges::begin(seq); it != std::ranges::end(seq); ++len, ++it)
        {
            // The rank cannot exceed 255 for single text and 254 for text collections as they are reserved as sentinels
            // for the indexed text.
            assert(seqan3::to_rank(static_cast<index_alphabet_type>(*it)) <
                   ((index_type::text_layout_mode == text_layout::single) ? 255 : 254));

            c = seqan3::to_rank(static_cast<index_alphabet_type>(*it)) + 1;

            new_parent_lb = _lb;
            new_parent_rb = _rb;
            if (!backward_search(index->index, c, _lb, _rb))
                return false;
        }

        parent_lb = new_parent_lb;
        parent_rb = new_parent_rb;
        node = {_lb, _rb, len + node.depth, c};
        return true;
    }

    /*!\brief Tries to replace the rightmost character of the query by the next lexicographically larger character such
     *        that the query is found in the text.
     *        \if DEV
     *            Moves the cursor to the right sibling of the current suffix tree node. It would be equivalent to
     *            going up an edge and going down that edge with the smallest character that is larger than the
     *            previous searched character. Calling cycle_back() on an cursor pointing to the root node is
     *            undefined behaviour!
     *        \endif
     * \returns `true` if there exists a query in the text where the rightmost character of the query is
     *          lexicographically larger than the current rightmost character of the query.
     *
     * Example:
     *
     * \include test/snippet/search/fm_index_cursor.cpp
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
        assert(index != nullptr && query_length() > 0);
        // parent_lb > parent_rb --> invalid interval
        assert(parent_lb <= parent_rb);

        sdsl_char_type c = node.last_char + 1;
        size_type _lb = parent_lb, _rb = parent_rb;

        while (c < sigma && !backward_search(index->index, index->index.comp2char[c], _lb, _rb))
        {
            ++c;
        }

        if (c != sigma) // Collection has additional sentinel as delimiter
        {
            node = {_lb, _rb, node.depth, c};
            return true;
        }
        return false;
    }

    /*!\brief Outputs the rightmost rank.
     * \returns Rightmost rank.
     *
     * Example:
     *
     * \include test/snippet/search/fm_index_cursor.cpp
     *
     * ### Complexity
     *
     * Constant.
     *
     * ### Exceptions
     *
     * No-throw guarantee.
     */
    size_type last_rank() const noexcept
    {
        // parent_lb > parent_rb --> invalid interval
        assert(index != nullptr && query_length() > 0 && parent_lb <= parent_rb);

        return index->index.comp2char[node.last_char] - 1; // text is not allowed to contain ranks of 0
    }

    /*!\brief Returns the length of the searched query.
     *        \if DEV
     *            Returns the depth of the cursor node in the implicit suffix tree.
     *        \endif
     * \returns Length of query.
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
        assert(node.depth != 0 || (node.lb == 0 && node.rb == index->size() - 1)); // depth == 0 -> root node

        return node.depth;
    }

    /*!\brief Returns the searched query.
     * \tparam text_t The type of the text used to build the index; must model std::ranges::input_range.
     * \param[in] text Text that was used to build the index.
     * \returns Searched query.
     *
     * \if DEV
     * Returns the concatenation of all edges from the root node to the cursors current node.
     * \endif
     *
     * ### Complexity
     *
     * \f$O(SAMPLING\_RATE * T_{BACKWARD\_SEARCH}) + query\_length()\f$
     *
     * ### Exceptions
     *
     * No-throw guarantee.
     */
    template <std::ranges::range text_t>
    auto path_label(text_t && text) const noexcept
    //!\cond
        requires index_t::text_layout_mode == text_layout::single
    //!\endcond
    {
        static_assert(std::ranges::input_range<text_t>, "The text must model input_range.");
        static_assert(dimension_v<text_t> == 1, "The input cannot be a text collection.");
        static_assert(std::same_as<innermost_value_type_t<text_t>, index_alphabet_type>,
                      "The alphabet types of the given text and index differ.");
        assert(index != nullptr);

        size_type const query_begin = offset() - index->index[node.lb];
        return text | views::slice(query_begin, query_begin + query_length());
    }

    //!\overload
    template <std::ranges::range text_t>
    auto path_label(text_t && text) const noexcept
    //!\cond
        requires index_t::text_layout_mode == text_layout::collection
    //!\endcond
    {
        static_assert(std::ranges::input_range<text_t>, "The text collection must model input_range.");
        static_assert(dimension_v<text_t> == 2, "The input must be a text collection.");
        static_assert(std::same_as<innermost_value_type_t<text_t>, index_alphabet_type>,
                      "The alphabet types of the given text and index differ.");
        assert(index != nullptr);

        // Position of query in concatenated text.
        size_type const location = offset() - index->index[node.lb];

        // The rank represents the number of start positions of the individual sequences/texts in the collection
        // before position `location + 1` and thereby also the number of delimiters.
        size_type const rank = index->text_begin_rs.rank(location + 1);
        assert(rank > 0);
        size_type const text_id = rank - 1;

        // The start location of the `text_id`-th text in the sequence (position of the `rank`-th 1 in the bitvector).
        size_type const start_location = index->text_begin_ss.select(rank);
        // Substract lengths of previous sequences.
        size_type const query_begin = location - start_location;

        // Take subtext, slice query out of it
        return text[text_id] | views::slice(query_begin, query_begin + query_length());
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
        assert(index != nullptr);

        return 1 + node.rb - node.lb;
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
        requires index_t::text_layout_mode == text_layout::single
    //!\endcond
    {
        assert(index != nullptr);

        std::vector<size_type> occ(count());
        for (size_type i = 0; i < occ.size(); ++i)
        {
            occ[i] = offset() - index->index[node.lb + i];
        }
        return occ;
    }

    //!\overload
    std::vector<std::pair<size_type, size_type>> locate() const
    //!\cond
        requires index_t::text_layout_mode == text_layout::collection
    //!\endcond
    {
        assert(index != nullptr);

        std::vector<std::pair<size_type, size_type>> occ;
        occ.reserve(count());
        for (size_type i = 0; i < count(); ++i)
        {
            size_type loc = offset() - index->index[node.lb + i];
            size_type sequence_rank = index->text_begin_rs.rank(loc + 1);
            size_type sequence_position = loc - index->text_begin_ss.select(sequence_rank);
            occ.emplace_back(sequence_rank - 1, sequence_position);
        }
        return occ;
    }

    /*!\brief Locates the occurrences of the searched query in the text on demand, i.e. a ranges::view is returned and
     *        every position is located once it is accessed.
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
        requires index_t::text_layout_mode == text_layout::single
    //!\endcond
    {
        assert(index != nullptr);

        return std::views::iota(node.lb, node.lb + count())
               | std::views::transform([*this, _offset = offset()] (auto sa_pos) { return _offset - index->index[sa_pos]; });
    }

    //!\overload
    auto lazy_locate() const
    //!\cond
        requires index_t::text_layout_mode == text_layout::collection
    //!\endcond
    {
        assert(index != nullptr);

        return std::views::iota(node.lb, node.lb + count())
               | std::views::transform([*this, _offset = offset()] (auto sa_pos) { return _offset - index->index[sa_pos]; })
               | std::views::transform([*this] (auto loc)
               {
                   size_type sequence_rank = index->text_begin_rs.rank(loc + 1);
                   size_type sequence_position = loc - index->text_begin_ss.select(sequence_rank);
                   return std::make_pair(sequence_rank-1, sequence_position);
               });
    }

};

//!\}

} // namespace seqan3
