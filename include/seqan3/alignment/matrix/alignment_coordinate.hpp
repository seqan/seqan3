// -----------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2019, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2019, Knut Reinert & MPI für molekulare Genetik
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at: https://github.com/seqan/seqan3/blob/master/LICENSE
// -----------------------------------------------------------------------------------------------------

/*!\file
 * \author Marcel Ehrhardt <marcel.ehrhardt AT fu-berlin.de>
 * \brief Provides seqan3::detail::alignment_coordinate.
 */

#pragma once

#include <type_traits>

#include <seqan3/core/detail/strong_type.hpp>
#include <seqan3/io/stream/debug_stream.hpp>
#include <seqan3/std/concepts>
#include <seqan3/std/iterator>

namespace seqan3::detail
{
/*!\brief A strong type for designated initialisation of the column index of the seqan3::detail::alignment_coordinate.
 * \ingroup alignment_matrix
 * \tparam index_t The type of the index; must model std::UnsignedIntegral.
 */
template <std::UnsignedIntegral index_t>
struct column_index_type : detail::strong_type<index_t, column_index_type<index_t>>
{
    //!!\brief Introduce base class constructor into this type's definition.
    using detail::strong_type<index_t, column_index_type<index_t>>::strong_type;
};

/*!\name Type deduction guides
 * \{
 */
//!\brief The default constructed element deduces to `size_t`.
column_index_type() -> column_index_type<size_t>;
//!\brief Deduces the position type from the constructor argument.
template <std::UnsignedIntegral index_t>
column_index_type(index_t const) -> column_index_type<size_t>;
//!\}

/*!\brief A strong type for designated initialisation of the row index of the seqan3::detail::alignment_coordinate.
 * \ingroup alignment_matrix
 * \tparam index_t The type of the index; must model std::UnsignedIntegral.
 */
template <std::UnsignedIntegral index_t>
struct row_index_type : detail::strong_type<index_t, row_index_type<index_t>>
{
    //!!\brief Introduce base class constructor into this type's definition.
    using detail::strong_type<index_t, row_index_type<index_t>>::strong_type;
};

/*!\name Type deduction guides
 * \{
 */
//!\brief The default constructed element deduces to `size_t`.
row_index_type() -> row_index_type<size_t>;
//!\brief Deduces the position type from the constructor argument.
template <std::UnsignedIntegral index_t>
row_index_type(index_t const) -> row_index_type<size_t>;
//!\}

/*!\brief Represents a state to specify the implementation of the seqan3::detail::advanceable_alignment_coordinate
 * \ingroup alignment_matrix
 *
 * \details
 *
 * The class seqan3::detail::alignment_coordinate can be extended with an incrementable and decrementable policy,
 * such that it can be used as a value type inside of a iota_view. This state offers three policies: none, which
 * leaves the functionality of seqan3::detail::alignment_coordinate untouched; column, which adds the respective
 * functionality only for the column index and row, which adds the respective functionality only for the row index.
 */
enum struct advanceable_alignment_coordinate_state : uint8_t
{
    //!\brief The corresponding alignment coordinate will not be incrementable/decrementable.
    none,
    //!\brief The corresponding alignment coordinate will be incrementable/decrementable in the column index.
    column,
    //!\brief The corresponding alignment coordinate will be incrementable/decrementable in the row index.
    row
};

/*!\brief Implements an internal alignment coordinate that can be used as an argument to the std::ranges::iota_view.
 * \ingroup alignment_matrix
 * \implements std::EqualityComparable
 * \tparam index_t  The type of the sequence index; must model std::UnsignedIntegral.
 * \tparam states   A non-type template flag to enable a specific advanceable policy.
 *                  See seqan3::detail::advanceable_alignment_coordinate_state
 *
 * \details
 *
 * This class provides all members to make the `advanceable_alignment_coordinate` be usable in a std::ranges::iota_view.
 * For the purpose of alignments modelling only incrementable and decrementable would fully suffice.
 * Unfortunately, the current range implementation does not preserve std::ranges::BidirectionalRange properties,
 * so we need to model the full advanceable concept in order to preserve the std::ranges::RandomAccessRange properties.
 * This however can be relaxed if the range implementation fully complies with the current standard draft for Ranges,
 * and increment and decrement would be enough.
 */
template <std::UnsignedIntegral index_t,
          advanceable_alignment_coordinate_state state = advanceable_alignment_coordinate_state::none>
class advanceable_alignment_coordinate
{
public:

    //!\name Member types
    //!\{
    //!\brief Defines the difference type to model the std::WeaklyIncrementable concept.
    using difference_type = std::make_signed_t<index_t>;
    //!\}

    /*!\name Constructor, destructor and assignment
     * \{
     */
    constexpr advanceable_alignment_coordinate() noexcept = default;
    constexpr advanceable_alignment_coordinate(advanceable_alignment_coordinate const &) noexcept = default;
    constexpr advanceable_alignment_coordinate(advanceable_alignment_coordinate &&) noexcept = default;
    constexpr advanceable_alignment_coordinate & operator=(advanceable_alignment_coordinate const &) noexcept = default;
    constexpr advanceable_alignment_coordinate & operator=(advanceable_alignment_coordinate &&) noexcept = default;
    ~advanceable_alignment_coordinate() noexcept = default;

    //!\brief Copy-constructs from another advanceable_alignment_coordinate with a different policy.
    template <typename _index_t, advanceable_alignment_coordinate_state _state>
        requires std::Same<_index_t, index_t> && !std::Same<_state, state>
    constexpr advanceable_alignment_coordinate(advanceable_alignment_coordinate<_index_t, _state> const & other) :
        first_seq_pos{other.first_seq_pos},
        second_seq_pos{other.second_seq_pos}
    {}

    //!\brief Move-constructs from another advanceable_alignment_coordinate with a different policy.
    template <typename _index_t, advanceable_alignment_coordinate_state _state>
        requires std::Same<_index_t, index_t> && !std::Same<_state, state>
    constexpr advanceable_alignment_coordinate(advanceable_alignment_coordinate<_index_t, _state> && other) :
        first_seq_pos{std::move(other.first_seq_pos)},
        second_seq_pos{std::move(other.second_seq_pos)}
    {}

    /*!\brief Save construction from column and row indices.
     * \param c_idx The respective column index within the matrix. Of type seqan3::detail::column_index_type.
     * \param r_idx The respective row index within the matrix. Of type seqan3::detail::row_index_type.
     */
    constexpr advanceable_alignment_coordinate(column_index_type<index_t> const c_idx,
                                               row_index_type<index_t> const r_idx) noexcept :
        first_seq_pos{c_idx.get()},
        second_seq_pos{r_idx.get()}
    {}
    //!\}

     //!\cond
     constexpr friend bool operator==(advanceable_alignment_coordinate const & lhs,
                                      advanceable_alignment_coordinate const & rhs) noexcept
     {
         return (lhs.first_seq_pos == rhs.first_seq_pos) &&
                (lhs.second_seq_pos == rhs.second_seq_pos);
     }

     constexpr friend bool operator!=(advanceable_alignment_coordinate const & lhs,
                                      advanceable_alignment_coordinate const & rhs) noexcept
     {
         return !(lhs == rhs);
     }
     //!\endcond

    /*!\name Member functions
     * \brief Advances or decrements the respective column or row coordinate depending on the set policy.
     *        These member functions are not available if the non-type template parameter `state` was set to
     *        seqan3::detail::advanceable_alignment_coordinate_state::none.
     * \{
     */

    constexpr advanceable_alignment_coordinate & operator++(/*pre-increment*/) noexcept
    //!\cond
        requires state != advanceable_alignment_coordinate_state::none
    //!\endcond
    {
        if constexpr (state == advanceable_alignment_coordinate_state::column)
            ++this->first_seq_pos;
        else
            ++this->second_seq_pos;
        return *this;
    }

    constexpr advanceable_alignment_coordinate operator++(int /*post-increment*/) noexcept
    //!\cond
        requires state != advanceable_alignment_coordinate_state::none
    //!\endcond
    {
        advanceable_alignment_coordinate tmp{*this};
        ++(*this);
        return tmp;
    }

    constexpr advanceable_alignment_coordinate & operator--(/*pre-decrement*/) noexcept
    //!\cond
        requires state != advanceable_alignment_coordinate_state::none
    //!\endcond
    {
        if constexpr (state == advanceable_alignment_coordinate_state::column)
            --this->first_seq_pos;
        else
            --this->second_seq_pos;
        return *this;
    }

    constexpr advanceable_alignment_coordinate operator--(int /*post-decrement*/) noexcept
    //!\cond
        requires state != advanceable_alignment_coordinate_state::none
    //!\endcond
    {
        advanceable_alignment_coordinate tmp{*this};
        --(*this);
        return tmp;
    }

    constexpr advanceable_alignment_coordinate & operator+=(difference_type const offset) noexcept
    //!\cond
        requires state != advanceable_alignment_coordinate_state::none
    //!\endcond
    {
        if constexpr (state == advanceable_alignment_coordinate_state::column)
            this->first_seq_pos += offset;
        else
            this->second_seq_pos += offset;
        return *this;
    }

    constexpr advanceable_alignment_coordinate & operator-=(difference_type const offset) noexcept
    //!\cond
        requires state != advanceable_alignment_coordinate_state::none
    //!\endcond
    {
        if constexpr (state == advanceable_alignment_coordinate_state::column)
            this->first_seq_pos -= offset;
        else
            this->second_seq_pos -= offset;
        return *this;
    }

    constexpr advanceable_alignment_coordinate operator+(difference_type const offset) const noexcept
    //!\cond
        requires state != advanceable_alignment_coordinate_state::none
    //!\endcond
    {
        advanceable_alignment_coordinate tmp{*this};
        tmp += offset;
        return tmp;
    }

    constexpr advanceable_alignment_coordinate operator-(difference_type const offset) const noexcept
    //!\cond
        requires state != advanceable_alignment_coordinate_state::none
    //!\endcond
    {
        advanceable_alignment_coordinate tmp{*this};
        tmp -= offset;
        return tmp;
    }

    constexpr difference_type operator-(advanceable_alignment_coordinate const & other) const noexcept
    //!\cond
        requires state != advanceable_alignment_coordinate_state::none
    //!\endcond
    {
        if constexpr (state == advanceable_alignment_coordinate_state::column)
            return this->first_seq_pos - other.first_seq_pos;
        else
            return this->second_seq_pos - other.second_seq_pos;
    }
    //!\}

    //!\brief Non-member function
    //!\{

    /*!\brief Advances the respective coordinate depending on the set policy by the given offset. This function is
              not available if the policy was set to seqan3::detail::advanceable_alignment_coordinate_state::none.
     */
    constexpr friend advanceable_alignment_coordinate operator+(difference_type const offset,
                                         advanceable_alignment_coordinate const & me) noexcept
    //!\cond
        requires state != advanceable_alignment_coordinate_state::none
    //!\endcond
    {
        return me + offset;
    }
    //!\}

    //!\brief The begin/end position of the alignment in the first sequence.
    index_t first_seq_pos{};
    //!\brief The begin/end position of the alignment in the second sequence.
    index_t second_seq_pos{};
};

/*!\name Type deduction guides
 * \{
 */
advanceable_alignment_coordinate() -> advanceable_alignment_coordinate<size_t>;
//!\brief Deduces the position type from the constructor.
template <std::UnsignedIntegral index_t>
advanceable_alignment_coordinate(column_index_type<index_t> const, row_index_type<index_t> const) ->
    advanceable_alignment_coordinate<size_t>;
//!\}

} // namespace seqan3::detail

namespace seqan3
{

/*!\brief Represents the begin/end of the pairwise alignment in the respective sequences.
 * \ingroup alignment_matrix
 *
 * \details
 * \if DEV
 * This class only gives access to the respective positions of the sequences and is meant for
 * the user interface. For the implementation of the alignment algorithms the
 * seqan3::detail::advanceable_alignment_coordinate is used.
 * \endif
 */
class alignment_coordinate
//!\cond DEV
    : public detail::advanceable_alignment_coordinate<size_t>
//!\endcond
{
    //!\cond DEV
    //!\brief The type of the base class.
    using base_t = detail::advanceable_alignment_coordinate<size_t>;
    //!\endcond

public:

    /*!\name Constructor, destructor and assignment
     * \{
     */
    constexpr alignment_coordinate() = default;
    constexpr alignment_coordinate(alignment_coordinate const &) = default;
    constexpr alignment_coordinate(alignment_coordinate &&) = default;
    constexpr alignment_coordinate & operator=(alignment_coordinate const &) = default;
    constexpr alignment_coordinate & operator=(alignment_coordinate &&) = default;
    ~alignment_coordinate() = default;

    //!\cond DEV
    //!\brief Constructs from the seqan3::detail::advanceable_alignment_coordinate base class.
    constexpr alignment_coordinate(base_t const & base) : base_t{base}
    {}

    //!\brief Constructs from the seqan3::detail::advanceable_alignment_coordinate base class.
    constexpr alignment_coordinate(base_t && base) : base_t{std::move(base)}
    {}

    /*!\brief Save construction from column and row indices.
     * \tparam index_t The index type used to store the column and row position.
     * \param[in] c_idx The respective column index within the matrix. Of type seqan3::detail::column_index_type.
     * \param[in] r_idx The respective row index within the matrix. Of type seqan3::detail::row_index_type.
     */
    template <typename index_t>
    constexpr alignment_coordinate(detail::column_index_type<index_t> const c_idx,
                                   detail::row_index_type<index_t> const r_idx) noexcept :
        base_t{c_idx, r_idx}
    {}
    //!\endcond
    //!\}

    //!\brief The begin/end position of the alignment in the first sequence.
    using base_t::first_seq_pos;
    //!\brief The begin/end position of the alignment in the second sequence.
    using base_t::second_seq_pos;

    //!\brief The begin/end position of the alignment in the first sequence.
    SEQAN3_DOXYGEN_ONLY(size_t first_seq_pos;)
    //!\brief The begin/end position of the alignment in the second sequence.
    SEQAN3_DOXYGEN_ONLY(size_t second_seq_pos;)
};

/*!\brief A seqan3::alignment_coordinate can be printed to the seqan3::debug_stream.
 * \tparam    coordinate_type The alignment coordinate type.
 * \param[in] s               The seqan3::debug_stream.
 * \param[in] c               The alignment coordinate to print.
 * \relates seqan3::debug_stream_type
 *
 * \details
 *
 * Prints the alignment coordinate as a tuple.
 */
template <typename coordinate_type>
//!\cond
    requires std::Same<remove_cvref_t<coordinate_type>, alignment_coordinate>
//!\endcond
inline debug_stream_type & operator<<(debug_stream_type & s, coordinate_type && c)
{
    s << std::tie(c.first_seq_pos, c.second_seq_pos);
    return s;
}

} // namespace seqan3
