// SPDX-FileCopyrightText: 2006-2025 Knut Reinert & Freie Universität Berlin
// SPDX-FileCopyrightText: 2016-2025 Knut Reinert & MPI für molekulare Genetik
// SPDX-License-Identifier: BSD-3-Clause

/*!\file
 * \author Marcel Ehrhardt <marcel.ehrhardt AT fu-berlin.de>
 * \brief Provides seqan3::detail::advanceable_alignment_coordinate.
 */

#pragma once

#include <concepts>
#include <iterator>
#include <type_traits>

#include <seqan3/alignment/matrix/detail/matrix_coordinate.hpp>
#include <seqan3/core/debug_stream/debug_stream_type.hpp>
#include <seqan3/core/debug_stream/tuple.hpp>
#include <seqan3/core/detail/strong_type.hpp>
#include <seqan3/core/detail/template_inspection.hpp>

namespace seqan3::detail
{

/*!\brief Represents a state to specify the implementation of the seqan3::detail::advanceable_alignment_coordinate.
 * \ingroup alignment_matrix
 *
 * \details
 *
 * The class seqan3::detail::advanceable_alignment_coordinate can be extended with an incrementable and decrementable
 * policy such that it can be used as a value type inside of a iota_view. This state offers three policies: none, which
 * leaves the functionality of seqan3::detail::advanceable_alignment_coordinate untouched; column, which adds the
 * respective functionality only for the column index and row, which adds the respective functionality only for the row
 * index.
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
    constexpr advanceable_alignment_coordinate() noexcept = default; //!< Defaulted
    //!\brief Defaulted
    constexpr advanceable_alignment_coordinate(advanceable_alignment_coordinate const &) noexcept = default;
    constexpr advanceable_alignment_coordinate(advanceable_alignment_coordinate &&) noexcept = default; //!< Defaulted
    //!\brief Defaulted
    constexpr advanceable_alignment_coordinate & operator=(advanceable_alignment_coordinate const &) noexcept = default;
    //!\brief Defaulted
    constexpr advanceable_alignment_coordinate & operator=(advanceable_alignment_coordinate &&) noexcept = default;
    ~advanceable_alignment_coordinate() noexcept = default; //!< Defaulted

    //!\brief Copy-constructs from another advanceable_alignment_coordinate with a different policy.
    template <advanceable_alignment_coordinate_state other_state>
        requires (other_state != state)
    constexpr advanceable_alignment_coordinate(advanceable_alignment_coordinate<other_state> const & other) :
        first{other.first},
        second{other.second}
    {}

    //!\brief Move-constructs from another advanceable_alignment_coordinate with a different policy.
    template <advanceable_alignment_coordinate_state other_state>
        requires (other_state != state)
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

    constexpr friend bool operator<(advanceable_alignment_coordinate const & lhs,
                                    advanceable_alignment_coordinate const & rhs) noexcept
    {
        return std::tie(lhs.first, lhs.second) < std::tie(rhs.first, rhs.second);
    }

    constexpr friend bool operator>=(advanceable_alignment_coordinate const & lhs,
                                     advanceable_alignment_coordinate const & rhs) noexcept
    {
        return std::tie(lhs.first, lhs.second) >= std::tie(rhs.first, rhs.second);
    }

    constexpr friend bool operator>(advanceable_alignment_coordinate const & lhs,
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
        requires (state != advanceable_alignment_coordinate_state::none)
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
        requires (state != advanceable_alignment_coordinate_state::none)
    {
        advanceable_alignment_coordinate tmp{*this};
        ++(*this);
        return tmp;
    }

    /*!\brief Decrements the coordinate depending on the set policy by one.
     * \return `*this`
     */
    constexpr advanceable_alignment_coordinate & operator--(/*pre-decrement*/) noexcept
        requires (state != advanceable_alignment_coordinate_state::none)
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
        requires (state != advanceable_alignment_coordinate_state::none)
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
        requires (state != advanceable_alignment_coordinate_state::none)
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
        requires (state != advanceable_alignment_coordinate_state::none)
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
        requires (state != advanceable_alignment_coordinate_state::none)
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
        requires (state != advanceable_alignment_coordinate_state::none)
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
        requires (state != advanceable_alignment_coordinate_state::none)
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
        requires (state != advanceable_alignment_coordinate_state::none)
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

/*!\brief The printer for seqan3::detail::advanceable_alignment_coordinate.
 *
 * Prints the alignment coordinate as a tuple of the column and row index.
 *
 * \tparam state_t The state of the detail::advanceable_alignment_coordinate.
 * \ingroup alignment_matrix
 */
template <auto state_t>
struct advanceable_alignment_coordinate_printer<detail::advanceable_alignment_coordinate<state_t>>
{
    /*!\brief The function call operator that prints the coordinate to the given stream.
     * \tparam stream_t The type of the stream.
     * \param[in,out] stream The stream to print to.
     * \param[in] arg The alignment coordinate to print.
     */
    template <typename stream_t>
    constexpr void operator()(stream_t & stream, detail::advanceable_alignment_coordinate<state_t> const arg) const
    {
        stream << std::tie(arg.first, arg.second);
    }
};

} // namespace seqan3
