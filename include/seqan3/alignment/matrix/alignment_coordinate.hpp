// -----------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2020, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2020, Knut Reinert & MPI für molekulare Genetik
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at: https://github.com/seqan/seqan3/blob/master/LICENSE.md
// -----------------------------------------------------------------------------------------------------

/*!\file
 * \author Marcel Ehrhardt <marcel.ehrhardt AT fu-berlin.de>
 * \brief Provides seqan3::detail::alignment_coordinate.
 */

#pragma once

#include <type_traits>

#include <seqan3/alignment/matrix/detail/matrix_coordinate.hpp>
#include <seqan3/core/detail/debug_stream_tuple.hpp>
#include <seqan3/core/detail/debug_stream_type.hpp>
#include <seqan3/core/detail/strong_type.hpp>
#include <seqan3/core/type_traits/template_inspection.hpp>
#include <seqan3/std/concepts>
#include <seqan3/std/iterator>

namespace seqan3::detail
{

/*!\brief Represents a state to specify the implementation of the seqan3::detail::advanceable_alignment_coordinate.
 * \ingroup alignment_matrix
 *
 * \details
 *
 * The class seqan3::detail::advanceable_alignment_coordinate can be extended with an incrementable and decrementable
 * policy such that it can be used as a value type inside of a iota_view. This state offers three policies: none, which
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
 * \implements std::equality_comparable
 * \tparam states   A non-type template flag to enable a specific advanceable policy.
 * \see seqan3::detail::advanceable_alignment_coordinate_state.
 *
 * \details
 *
 * This class provides all members to make the `advanceable_alignment_coordinate` be usable in a std::ranges::iota_view.
 * For the purpose of alignments, modelling only incrementable and decrementable would fully suffice.
 * Unfortunately, the current range implementation does not preserve std::ranges::bidirectional_range properties,
 * so we need to model the full advanceable concept in order to preserve the std::ranges::random_access_range properties.
 * This, however, can be relaxed if the range implementation fully complies with the current standard draft for Ranges,
 * as increment and decrement would be enough and std::views::iota would preserve std::ranges::bidirectional_range.
 */
template <advanceable_alignment_coordinate_state state = advanceable_alignment_coordinate_state::none>
class advanceable_alignment_coordinate
{
public:

    /*!\name Member types
     * \{
     */

    //!\brief Defines the difference type to model the std::weakly_incrementable concept.
    using difference_type = std::make_signed_t<size_t>;
    //!\}

    /*!\name Constructors, destructor and assignment
     * \{
     */
    constexpr advanceable_alignment_coordinate() noexcept = default;                                    //!< Defaulted
    //!\brief Defaulted
    constexpr advanceable_alignment_coordinate(advanceable_alignment_coordinate const &) noexcept = default;
    constexpr advanceable_alignment_coordinate(advanceable_alignment_coordinate &&) noexcept = default; //!< Defaulted
    //!\brief Defaulted
    constexpr advanceable_alignment_coordinate & operator=(advanceable_alignment_coordinate const &) noexcept = default;
    //!\brief Defaulted
    constexpr advanceable_alignment_coordinate & operator=(advanceable_alignment_coordinate &&) noexcept = default;
    ~advanceable_alignment_coordinate() noexcept = default;                                             //!< Defaulted

    //!\brief Copy-constructs from another advanceable_alignment_coordinate with a different policy.
    template <advanceable_alignment_coordinate_state other_state>
    //!\cond
        requires other_state != state
    //!\endcond
    constexpr advanceable_alignment_coordinate(advanceable_alignment_coordinate<other_state> const & other) :
        first{other.first},
        second{other.second}
    {}

    //!\brief Move-constructs from another advanceable_alignment_coordinate with a different policy.
    template <advanceable_alignment_coordinate_state other_state>
    //!\cond
        requires other_state != state
    //!\endcond
    constexpr advanceable_alignment_coordinate(advanceable_alignment_coordinate<other_state> && other) :
        first{std::move(other.first)},
        second{std::move(other.second)}
    {}

    /*!\brief Construction from the respective column and row indices.
     * \param c_idx The respective column index within the matrix. Of type seqan3::detail::column_index_type.
     * \param r_idx The respective row index within the matrix. Of type seqan3::detail::row_index_type.
     */
    constexpr advanceable_alignment_coordinate(column_index_type<size_t> const c_idx,
                                               row_index_type<size_t> const r_idx) noexcept :
        first{c_idx.get()},
        second{r_idx.get()}
    {}
    //!\}

     //!\cond
    constexpr friend bool operator==(advanceable_alignment_coordinate const & lhs,
                                     advanceable_alignment_coordinate const & rhs) noexcept
    {
        return std::tie(lhs.first, lhs.second) == std::tie(rhs.first, rhs.second);
    }

    constexpr friend bool operator!=(advanceable_alignment_coordinate const & lhs,
                                     advanceable_alignment_coordinate const & rhs) noexcept
    {
        return std::tie(lhs.first, lhs.second) != std::tie(rhs.first, rhs.second);
    }

    constexpr friend bool operator<=(advanceable_alignment_coordinate const & lhs,
                                     advanceable_alignment_coordinate const & rhs) noexcept
    {
        return std::tie(lhs.first, lhs.second) <= std::tie(rhs.first, rhs.second);
    }

    constexpr friend bool operator< (advanceable_alignment_coordinate const & lhs,
                                     advanceable_alignment_coordinate const & rhs) noexcept
    {
        return std::tie(lhs.first, lhs.second) < std::tie(rhs.first, rhs.second);
    }

    constexpr friend bool operator>=(advanceable_alignment_coordinate const & lhs,
                                     advanceable_alignment_coordinate const & rhs) noexcept
    {
        return std::tie(lhs.first, lhs.second) >= std::tie(rhs.first, rhs.second);
    }

    constexpr friend bool operator> (advanceable_alignment_coordinate const & lhs,
                                     advanceable_alignment_coordinate const & rhs) noexcept
    {
        return std::tie(lhs.first, lhs.second) > std::tie(rhs.first, rhs.second);
    }
     //!\endcond

    /*!\name Member functions
     * \brief Advances or decrements the respective column or row coordinate depending on the set policy.
     *        These member functions are not available if the non-type template parameter `state` was set to
     *        seqan3::detail::advanceable_alignment_coordinate_state::none.
     * \{
     */

    /*!\brief Increments the coordinate depending on the set policy by one.
     * \return `*this`
     */
    constexpr advanceable_alignment_coordinate & operator++(/*pre-increment*/) noexcept
    //!\cond
        requires state != advanceable_alignment_coordinate_state::none
    //!\endcond
    {
        if constexpr (state == advanceable_alignment_coordinate_state::column)
            ++this->first;
        else
            ++this->second;
        return *this;
    }

    /*!\brief Post-increments the coordinate depending on the set policy by one.
     * \return a seqan3::detail::advanceable_alignment_coordinate that holds an unchanged value.
     */
    constexpr advanceable_alignment_coordinate operator++(int /*post-increment*/) noexcept
    //!\cond
        requires state != advanceable_alignment_coordinate_state::none
    //!\endcond
    {
        advanceable_alignment_coordinate tmp{*this};
        ++(*this);
        return tmp;
    }

    /*!\brief Decrements the coordinate depending on the set policy by one.
     * \return `*this`
     */
    constexpr advanceable_alignment_coordinate & operator--(/*pre-decrement*/) noexcept
    //!\cond
        requires state != advanceable_alignment_coordinate_state::none
    //!\endcond
    {
        if constexpr (state == advanceable_alignment_coordinate_state::column)
            --this->first;
        else
            --this->second;
        return *this;
    }

    /*!\brief Post-decrements the coordinate depending on the set policy by one.
     * \return a seqan3::detail::advanceable_alignment_coordinate that holds an unchanged value.
     */
    constexpr advanceable_alignment_coordinate operator--(int /*post-decrement*/) noexcept
    //!\cond
        requires state != advanceable_alignment_coordinate_state::none
    //!\endcond
    {
        advanceable_alignment_coordinate tmp{*this};
        --(*this);
        return tmp;
    }

    /*!\brief Returns the coordinate which is advanced depending on the set policy by `offset`.
     * \param offset The value to add to the coordinate.
     * \return `*this`
     */
    constexpr advanceable_alignment_coordinate & operator+=(difference_type const offset) noexcept
    //!\cond
        requires state != advanceable_alignment_coordinate_state::none
    //!\endcond
    {
        if constexpr (state == advanceable_alignment_coordinate_state::column)
            this->first += offset;
        else
            this->second += offset;
        return *this;
    }

    /*!\brief Returns the coordinate which is advanced depending on the set policy by -`offset`.
     * \param offset The value to subtract from the coordinate.
     * \return `*this`
     */
    constexpr advanceable_alignment_coordinate & operator-=(difference_type const offset) noexcept
    //!\cond
        requires state != advanceable_alignment_coordinate_state::none
    //!\endcond
    {
        if constexpr (state == advanceable_alignment_coordinate_state::column)
            this->first -= offset;
        else
            this->second -= offset;
        return *this;
    }

    /*!\brief Returns a new coordinate which is advanced depending on the set policy by `offset`.
     * \param offset The value to add to the coordinate.
     * \return a seqan3::detail::advanceable_alignment_coordinate that holds the updated value.
     */
    constexpr advanceable_alignment_coordinate operator+(difference_type const offset) const noexcept
    //!\cond
        requires state != advanceable_alignment_coordinate_state::none
    //!\endcond
    {
        advanceable_alignment_coordinate tmp{*this};
        tmp += offset;
        return tmp;
    }

    /*!\brief Returns a new coordinate which is advanced depending on the set policy by -`offset`.
     * \param offset The value to subtract from the coordinate.
     * \return a seqan3::detail::advanceable_alignment_coordinate that holds the updated value.
     */
    constexpr advanceable_alignment_coordinate operator-(difference_type const offset) const noexcept
    //!\cond
        requires state != advanceable_alignment_coordinate_state::none
    //!\endcond
    {
        advanceable_alignment_coordinate tmp{*this};
        tmp -= offset;
        return tmp;
    }

    /*!\brief Returns the difference of this and another coordinate depending on the set policy.
     * \param other The other coordinate.
     * \return the difference between the coordinates.
     */
    constexpr difference_type operator-(advanceable_alignment_coordinate const & other) const noexcept
    //!\cond
        requires state != advanceable_alignment_coordinate_state::none
    //!\endcond
    {
        if constexpr (state == advanceable_alignment_coordinate_state::column)
            return this->first - other.first;
        else
            return this->second - other.second;
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

    //!\brief The front/back position of the alignment in the first sequence.
    size_t first{};
    //!\brief The front/back position of the alignment in the second sequence.
    size_t second{};
};

} // namespace seqan3::detail

namespace seqan3
{

/*!\brief Represents the begin/end of the pairwise alignment in the respective sequences.
 * \ingroup alignment_matrix
 *
 * \if DEV
 * \details
 * This class only gives access to the respective positions of the sequences and is meant for
 * the user interface. The additional complexity of an advanceable coorindate using the
 * seqan3::detail::advanceable_alignment_coordinate is only necessary for the implementation of the pairwise
 * alignment algorithm. Within in the algorithm the coordinate is used in combination with a seqan3::views::iota to
 * keep track of the current position within the alignment matrix. For the user, however, this interface adds no
 * benefit as they are only interested in the front/back coordinates for the respective alignment.
 * \endif
 */
class alignment_coordinate
//!\cond DEV
    : public detail::advanceable_alignment_coordinate<detail::advanceable_alignment_coordinate_state::none>
//!\endcond
{
    //!\cond DEV
    //!\brief The type of the base class.
    using base_t = detail::advanceable_alignment_coordinate<detail::advanceable_alignment_coordinate_state::none>;
    //!\endcond

public:

    /*!\name Constructors, destructor and assignment
     * \{
     */
    constexpr alignment_coordinate() = default;                                         //!< Defaulted
    constexpr alignment_coordinate(alignment_coordinate const &) = default;             //!< Defaulted
    constexpr alignment_coordinate(alignment_coordinate &&) = default;                  //!< Defaulted
    constexpr alignment_coordinate & operator=(alignment_coordinate const &) = default; //!< Defaulted
    constexpr alignment_coordinate & operator=(alignment_coordinate &&) = default;      //!< Defaulted
    ~alignment_coordinate() = default;                                                  //!< Defaulted

    //!\cond DEV
    //!\brief Inherit the constructor from the base class.
    using base_t::base_t;

    //!\brief Constructs from the seqan3::detail::advanceable_alignment_coordinate base class.
    constexpr alignment_coordinate(base_t const & base) : base_t{base}
    {}

    //!\brief Constructs from the seqan3::detail::advanceable_alignment_coordinate base class.
    constexpr alignment_coordinate(base_t && base) : base_t{std::move(base)}
    {}
    //!\endcond
    //!\}

    using base_t::first;
    using base_t::second;

    //!\brief The begin/end position of the alignment in the first sequence.
    SEQAN3_DOXYGEN_ONLY(size_t first;)
    //!\brief The begin/end position of the alignment in the second sequence.
    SEQAN3_DOXYGEN_ONLY(size_t second;)

    //!\privatesection
    //!\brief Implicit conversion to seqan3::detail::matrix_coordinate.
    constexpr operator detail::matrix_coordinate() const
    {
        return detail::matrix_coordinate{detail::row_index_type{second}, detail::column_index_type{first}};
    }
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
template <typename coordinate_type, typename char_t>
//!\cond
    requires std::same_as<remove_cvref_t<coordinate_type>, alignment_coordinate> ||
             detail::is_value_specialisation_of_v<remove_cvref_t<coordinate_type>,
                                                  detail::advanceable_alignment_coordinate>
//!\endcond
inline debug_stream_type<char_t> & operator<<(debug_stream_type<char_t> & s, coordinate_type && c)
{
    s << std::tie(c.first, c.second);
    return s;
}

} // namespace seqan3
