// -----------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2020, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2020, Knut Reinert & MPI für molekulare Genetik
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at: https://github.com/seqan/seqan3/blob/master/LICENSE.md
// -----------------------------------------------------------------------------------------------------

/*!\file
 * \brief Provides seqan3::detail::affine_cell_proxy.
 * \author Rene Rahn <rene.rahn AT fu-berlin.de>
 */

#pragma once

#include <tuple>
#include <type_traits>

#include <seqan3/core/concept/core_language.hpp>
#include <seqan3/core/simd/concept.hpp>
#include <seqan3/core/tuple_utility.hpp>

namespace seqan3::detail
{

/*!\interface seqan3::detail::arithmetic_or_simd <>
 * \brief The concept for a type that models either seqan3::arithmetic or seqan3::simd::simd_concept.
 * \ingroup alignment_matrix
 */
//!\cond
template <typename t>
SEQAN3_CONCEPT arithmetic_or_simd = arithmetic<t> || simd_concept<t>;
//!\endcond

/*!\interface seqan3::detail::affine_score_cell <>
 * \extends seqan3::tuple_like
 * \brief The concept for a type that models an affine cell of the score matrix.
 * \ingroup alignment_matrix
 *
 * \details
 *
 * This concept describes the requirements an alignment matrix cell must fulfil to represent an affine score
 * matrix entry.
 */
//!\cond
template <typename t>
SEQAN3_CONCEPT affine_score_cell = tuple_like<t> &&
                                   std::tuple_size_v<t> == 3 &&
                                   arithmetic_or_simd<std::remove_reference_t<std::tuple_element_t<0, t>>> &&
                                   arithmetic_or_simd<std::remove_reference_t<std::tuple_element_t<1, t>>> &&
                                   arithmetic_or_simd<std::remove_reference_t<std::tuple_element_t<2, t>>>;
//!\endcond

/*!\brief A wrapper for an affine score matrix cell.
 * \implements seqan3::tuple_like
 * \ingroup alignment_matrix
 *
 * \tparam tuple_t The underlying cell type of the affine alignment matrix; must model
 *                 seqan3::detail::affine_score_cell.
 *
 * \details
 *
 * This wrapper provides a uniform access to the different elements of the cell within an affine score matrix. This
 * includes the optimal score, the horizontal gap score and the vertical gap score.
 */
template <affine_score_cell tuple_t>
//!\cond
    requires (std::tuple_size_v<tuple_t> == 3)
//!\endcond
class affine_cell_proxy : public tuple_t
{
public:
    /*!\name Constructors, destructor and assignment
     * \{
     */
    affine_cell_proxy()= default; //!< Defaulted.
    affine_cell_proxy(affine_cell_proxy const &) = default; //!< Defaulted.
    affine_cell_proxy(affine_cell_proxy &&) = default; //!< Defaulted.
    affine_cell_proxy & operator=(affine_cell_proxy const &) = default; //!< Defaulted.
    affine_cell_proxy & operator=(affine_cell_proxy &&) = default; //!< Defaulted.
    ~affine_cell_proxy() = default; //!< Defaulted.

    // Inherit the tuple constructors and assignment.
    using tuple_t::tuple_t;
    using tuple_t::operator=;
    //!\}

    /*!\name Score value accessor
     * \brief Specific accessor function to get the respective score value from an affine matrix cell.
     * \{
     */
    //!\brief Access the best score of the wrapped score matrix cell.
    decltype(auto) best_score() & noexcept { return get_score_impl<0>(*this); }
    //!\overload
    decltype(auto) best_score() const & noexcept { return get_score_impl<0>(*this); }
    //!\overload
    decltype(auto) best_score() && noexcept { return get_score_impl<0>(std::move(*this)); }
    //!\overload
    decltype(auto) best_score() const && noexcept
    { //Unfortunately gcc7 does not preserve the const && type: https://gcc.gnu.org/bugzilla/show_bug.cgi?id=94967
        using return_t = std::tuple_element_t<0, tuple_t>;
        return static_cast<return_t const &&>(get_score_impl<0>(std::move(*this)));
    }

    //!\brief Access the horizontal score of the wrapped score matrix cell.
    decltype(auto) horizontal_score() & noexcept { return get_score_impl<1>(*this); }
    //!\overload
    decltype(auto) horizontal_score() const & noexcept { return get_score_impl<1>(*this); }
    //!\overload
    decltype(auto) horizontal_score() && noexcept { return get_score_impl<1>(std::move(*this)); }
    //!\overload
     decltype(auto) horizontal_score() const && noexcept
    { //Unfortunately gcc7 does not preserve the const && type: https://gcc.gnu.org/bugzilla/show_bug.cgi?id=94967
        using return_t = std::tuple_element_t<1, tuple_t>;
        return static_cast<return_t const &&>(get_score_impl<1>(std::move(*this)));
    }

    //!\brief Access the vertical score of the wrapped score matrix cell.
    decltype(auto) vertical_score() & noexcept { return get_score_impl<2>(*this); }
    //!\overload
    decltype(auto) vertical_score() const & noexcept { return get_score_impl<2>(*this); }
    //!\overload
    decltype(auto) vertical_score() && noexcept { return get_score_impl<2>(std::move(*this)); }
    //!\overload
    decltype(auto) vertical_score() const && noexcept
    { //Unfortunately gcc7 does not preserve the const && type: https://gcc.gnu.org/bugzilla/show_bug.cgi?id=94967
        using return_t = std::tuple_element_t<2, tuple_t>;
        return static_cast<return_t const &&>(get_score_impl<2>(std::move(*this)));
    }
    //!\}

private:
    /*!\brief Implements the get interface for the various calls to receive the score value.
     * \tparam index The index of the tuple element to get; must be smaller than 3.
     * \tparam this_t The perfectly forwarded type of `*this`.
     *
     * \param[in] me The instance of `*this`.
     *
     * \returns The score value from the given tuple index.
     */
    template <size_t index, typename this_t>
    //!\cond
        requires (index < 3)
    //!\endcond
    static constexpr decltype(auto) get_score_impl(this_t && me) noexcept
    {
        using std::get;

        return get<index>(std::forward<this_t>(me));
    }
};
} // namespace seqan3::detail

namespace std
{
//!\cond
template <typename tuple_t>
struct tuple_size<seqan3::detail::affine_cell_proxy<tuple_t>> : tuple_size<tuple_t>
{};

template <size_t index, typename tuple_t>
struct tuple_element<index, seqan3::detail::affine_cell_proxy<tuple_t>> : tuple_element<index, tuple_t>
{};
//!\endcond
} // namespace std
