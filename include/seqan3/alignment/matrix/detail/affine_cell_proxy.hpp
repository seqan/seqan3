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

#include <seqan3/std/concepts>
#include <tuple>
#include <type_traits>

#include <seqan3/alignment/matrix/trace_directions.hpp>
#include <seqan3/core/detail/empty_type.hpp>
#include <seqan3/utility/detail/exposition_only_concept.hpp>
#include <seqan3/utility/simd/concept.hpp>
#include <seqan3/utility/tuple/concept.hpp>

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

/*!\interface seqan3::detail::tracedirections_or_simd <>
 * \brief The concept for a type that either is the same type as seqan3::detail::trace_directions or
 *        models the seqan3::simd::simd_concept.
 * \ingroup alignment_matrix
 */
//!\cond
template <typename t>
SEQAN3_CONCEPT tracedirections_or_simd = std::same_as<std::remove_cvref_t<t>, trace_directions> || simd_concept<t>;
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

/*!\interface seqan3::detail::affine_trace_cell <>
 * \extends seqan3::tuple_like
 * \brief The concept for a type that models an affine cell of the trace matrix.
 * \ingroup alignment_matrix
 *
 * \details
 *
 * This concept describes the requirements an alignment matrix cell must fulfil to represent an affine trace
 * matrix entry.
 */
//!\cond
template <typename t>
SEQAN3_CONCEPT affine_trace_cell = tuple_like<t> &&
                                   std::tuple_size_v<t> == 3 &&
                                   tracedirections_or_simd<std::remove_reference_t<std::tuple_element_t<0, t>>> &&
                                   tracedirections_or_simd<std::remove_reference_t<std::tuple_element_t<1, t>>> &&
                                   tracedirections_or_simd<std::remove_reference_t<std::tuple_element_t<2, t>>>;
//!\endcond

/*!\interface seqan3::detail::affine_score_and_trace_cell <>
 * \extends seqan3::tuple_like
 * \brief The concept for a type that models an affine cell of the combined score and trace matrix.
 * \ingroup alignment_matrix
 *
 * \details
 *
 * This concept describes the requirements an alignment matrix cell must fulfil to represent an affine score
 * matrix entry with the score and trace information.
 */
//!\cond
template <typename t>
SEQAN3_CONCEPT affine_score_and_trace_cell = tuple_like<t> &&
                                             std::tuple_size_v<t> == 2 &&
                                             affine_score_cell<std::tuple_element_t<0, t>> &&
                                             affine_trace_cell<std::tuple_element_t<1, t>>;
//!\endcond

/*!\brief A proxy for an affine score matrix cell.
 * \implements seqan3::tuple_like
 * \ingroup alignment_matrix
 *
 * \tparam tuple_t The underlying cell type of the affine alignment matrix; must model
 *                 seqan3::detail::affine_score_cell or seqan3::detail::affine_score_and_trace_cell.
 *
 * \details
 *
 * This wrapper provides a uniform access to the different elements of the cell within an affine score matrix. This
 * includes the best score, the horizontal gap score and the vertical gap score. In case of a combined alignment
 * matrix including the trace matrix, the interface is extended to also access the best, horizontal, and vertical trace
 * value.
 */
template <typename tuple_t>
//!\cond
    requires (affine_score_cell<tuple_t> || affine_score_and_trace_cell<tuple_t>)
//!\endcond
class affine_cell_proxy : public tuple_t
{
private:
    //!\brief The type of the score cell.
    using score_cell_type = std::conditional_t<affine_score_cell<tuple_t>, tuple_t, std::tuple_element_t<0, tuple_t>>;
    //!\brief The type of the trace cell (might be seqan3::detail::empty_type if not defined).
    using trace_cell_type = std::conditional_t<affine_score_and_trace_cell<tuple_t>,
                                               std::tuple_element_t<1, tuple_t>,
                                               empty_type>;

public:
    /*!\name Constructors, destructor and assignment
     * \{
     */
    affine_cell_proxy() = default; //!< Defaulted.
    affine_cell_proxy(affine_cell_proxy const &) = default; //!< Defaulted.
    affine_cell_proxy(affine_cell_proxy &&) = default; //!< Defaulted.
    affine_cell_proxy & operator=(affine_cell_proxy const &) = default; //!< Defaulted.
    affine_cell_proxy & operator=(affine_cell_proxy &&) = default; //!< Defaulted.
    ~affine_cell_proxy() = default; //!< Defaulted.

    // Inherit the base class's constructor to enable element-wise initialisation (direct and converting constructor).
    using tuple_t::tuple_t;

    //!\brief Converting constructor. Initialises from another tuple type.
    template <typename other_tuple_t>
    //!\cond
        requires std::constructible_from<tuple_t, other_tuple_t &&>
    //!\endcond
    explicit affine_cell_proxy(other_tuple_t && other) :
        tuple_t{std::forward<other_tuple_t>(other)}
    {}

    //!\brief Converting copy-constructor.
    template <typename other_tuple_t>
    //!\cond
        requires std::constructible_from<tuple_t, other_tuple_t const &>
    //!\endcond
    explicit affine_cell_proxy(affine_cell_proxy<other_tuple_t> const & other) :
        tuple_t{static_cast<other_tuple_t const &>(other)}
    {}

    //!\brief Converting move-constructor.
    template <typename other_tuple_t>
    //!\cond
        requires std::constructible_from<tuple_t, other_tuple_t>
    //!\endcond
    explicit affine_cell_proxy(affine_cell_proxy<other_tuple_t> && other) :
        tuple_t{static_cast<other_tuple_t &&>(std::move(other))}
    {}

    //!\brief Converting assignment. Initialises from another tuple type.
    template <typename other_tuple_t>
    //!\cond
        requires std::assignable_from<tuple_t &, other_tuple_t &&>
    //!\endcond
    affine_cell_proxy & operator=(other_tuple_t && other)
    {
        as_base() = std::forward<other_tuple_t>(other);
        return *this;
    }

    //!\brief Converting copy-assignment.
    template <typename other_tuple_t>
    //!\cond
        requires std::assignable_from<tuple_t &, other_tuple_t const &>
    //!\endcond
    affine_cell_proxy & operator=(affine_cell_proxy<other_tuple_t> const & other)
    {
        as_base() = static_cast<other_tuple_t const &>(other);
        return *this;
    }

    //!\brief Converting move-assignment.
    template <typename other_tuple_t>
    //!\cond
        requires std::assignable_from<tuple_t &, other_tuple_t>
    //!\endcond
    affine_cell_proxy & operator=(affine_cell_proxy<other_tuple_t> && other)
    {
        as_base() = static_cast<other_tuple_t &&>(std::move(other));
        return *this;
    }
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
    {
#if SEQAN3_WORKAROUND_GCC_94967
        // A simple std::move(...) does not work, because it would mess up tuple_element types like `int const &`
        using return_t = std::tuple_element_t<0, score_cell_type>;
        return static_cast<return_t const &&>(get_score_impl<0>(std::move(*this)));
#else // ^^^ workaround / no workaround vvv
        return get_score_impl<0>(std::move(*this));
#endif // SEQAN3_WORKAROUND_GCC_94967
    }

    //!\brief Access the horizontal score of the wrapped score matrix cell.
    decltype(auto) horizontal_score() & noexcept { return get_score_impl<1>(*this); }
    //!\overload
    decltype(auto) horizontal_score() const & noexcept { return get_score_impl<1>(*this); }
    //!\overload
    decltype(auto) horizontal_score() && noexcept { return get_score_impl<1>(std::move(*this)); }
    //!\overload
    decltype(auto) horizontal_score() const && noexcept
    {
#if SEQAN3_WORKAROUND_GCC_94967
        // A simple std::move(...) does not work, because it would mess up tuple_element types like `int const &`
        using return_t = std::tuple_element_t<1, score_cell_type>;
        return static_cast<return_t const &&>(get_score_impl<1>(std::move(*this)));
#else // ^^^ workaround / no workaround vvv
        return get_score_impl<1>(std::move(*this));
#endif // SEQAN3_WORKAROUND_GCC_94967
    }

    //!\brief Access the vertical score of the wrapped score matrix cell.
    decltype(auto) vertical_score() & noexcept { return get_score_impl<2>(*this); }
    //!\overload
    decltype(auto) vertical_score() const & noexcept { return get_score_impl<2>(*this); }
    //!\overload
    decltype(auto) vertical_score() && noexcept { return get_score_impl<2>(std::move(*this)); }
    //!\overload
    decltype(auto) vertical_score() const && noexcept
    {
#if SEQAN3_WORKAROUND_GCC_94967
        // A simple std::move(...) does not work, because it would mess up tuple_element types like `int const &`
        using return_t = std::tuple_element_t<2, score_cell_type>;
        return static_cast<return_t const &&>(get_score_impl<2>(std::move(*this)));
#else // ^^^ workaround / no workaround vvv
        return get_score_impl<2>(std::move(*this));
#endif // SEQAN3_WORKAROUND_GCC_94967
    }
    //!\}

    /*!\name Trace value accessor
     * \brief Specific accessor function to get the respective trace value from an affine matrix cell.
     * \{
     */
    //!\brief Access the optimal score of the wrapped score matrix cell.
    decltype(auto) best_trace() & noexcept
    //!\cond
        requires affine_score_and_trace_cell<tuple_t>
    //!\endcond
    {
        return get_trace_impl<0>(*this);
    }
    //!\overload
    decltype(auto) best_trace() const & noexcept
    //!\cond
        requires affine_score_and_trace_cell<tuple_t>
    //!\endcond
    {
        return get_trace_impl<0>(*this);
    }
    //!\overload
    decltype(auto) best_trace() && noexcept
    //!\cond
        requires affine_score_and_trace_cell<tuple_t>
    //!\endcond
    {
        return get_trace_impl<0>(std::move(*this));
    }
    //!\overload
    decltype(auto) best_trace() const && noexcept
    //!\cond
        requires affine_score_and_trace_cell<tuple_t>
    //!\endcond
    {
#if SEQAN3_WORKAROUND_GCC_94967
        // A simple std::move(...) does not work, because it would mess up tuple_element types like `int const &`
        using return_t = std::tuple_element_t<0, trace_cell_type>;
        return static_cast<return_t const &&>(get_trace_impl<0>(std::move(*this)));
#else // ^^^ workaround / no workaround vvv
        return get_trace_impl<0>(std::move(*this));
#endif // SEQAN3_WORKAROUND_GCC_94967
    }

    //!\brief Access the horizontal score of the wrapped score matrix cell.
    decltype(auto) horizontal_trace() & noexcept
    //!\cond
        requires affine_score_and_trace_cell<tuple_t>
    //!\endcond
    {
        return get_trace_impl<1>(*this);
    }
    //!\overload
    decltype(auto) horizontal_trace() const & noexcept
    //!\cond
        requires affine_score_and_trace_cell<tuple_t>
    //!\endcond
    {
        return get_trace_impl<1>(*this);
    }
    //!\overload
    decltype(auto) horizontal_trace() && noexcept
    //!\cond
        requires affine_score_and_trace_cell<tuple_t>
    //!\endcond
    {
        return get_trace_impl<1>(std::move(*this));
    }
    //!\overload
    decltype(auto) horizontal_trace() const && noexcept
    //!\cond
        requires affine_score_and_trace_cell<tuple_t>
    //!\endcond
    {
#if SEQAN3_WORKAROUND_GCC_94967
        // A simple std::move(...) does not work, because it would mess up tuple_element types like `int const &`
        using return_t = std::tuple_element_t<1, trace_cell_type>;
        return static_cast<return_t const &&>(get_trace_impl<1>(std::move(*this)));
#else // ^^^ workaround / no workaround vvv
        return get_trace_impl<1>(std::move(*this));
#endif // SEQAN3_WORKAROUND_GCC_94967
    }

    //!\brief Access the vertical score of the wrapped score matrix cell.
    decltype(auto) vertical_trace() & noexcept
    //!\cond
        requires affine_score_and_trace_cell<tuple_t>
    //!\endcond
    {
        return get_trace_impl<2>(*this);
    }
    //!\overload
    decltype(auto) vertical_trace() const & noexcept
    //!\cond
        requires affine_score_and_trace_cell<tuple_t>
    //!\endcond
    {
        return get_trace_impl<2>(*this);
    }
    //!\overload
    decltype(auto) vertical_trace() && noexcept
    //!\cond
        requires affine_score_and_trace_cell<tuple_t>
    //!\endcond
    {
        return get_trace_impl<2>(std::move(*this));
    }
    //!\overload
    decltype(auto) vertical_trace() const && noexcept
    //!\cond
        requires affine_score_and_trace_cell<tuple_t>
    //!\endcond
    {
#if SEQAN3_WORKAROUND_GCC_94967
        // A simple std::move(...) does not work, because it would mess up tuple_element types like `int const &`
        using return_t = std::tuple_element_t<2, trace_cell_type>;
        return static_cast<return_t const &&>(get_trace_impl<2>(std::move(*this)));
#else // ^^^ workaround / no workaround vvv
        return get_trace_impl<2>(std::move(*this));
#endif // SEQAN3_WORKAROUND_GCC_94967
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

        if constexpr (affine_score_cell<tuple_t>)
            return get<index>(std::forward<this_t>(me));
        else
            return get<index>(get<0>(std::forward<this_t>(me)));
    }

    /*!\brief Implements the get interface for the various calls to receive the trace value.
     * \tparam index The index of the tuple element to get; must be smaller than 3.
     * \tparam this_t The perfectly forwarded type of `*this`.
     *
     * \param[in] me The instance of `*this`.
     *
     * \returns The trace value from the given tuple index.
     */
    template <size_t index, typename this_t>
    //!\cond
        requires (index < 3 && affine_score_and_trace_cell<tuple_t>)
    //!\endcond
    static constexpr decltype(auto) get_trace_impl(this_t && me) noexcept
    {
        using std::get;

        return get<index>(get<1>(std::forward<this_t>(me)));
    }

    //!\brief Casts `this` to the base class type.
    tuple_t & as_base() & noexcept
    {
        return static_cast<tuple_t &>(*this);
    }
};
} // namespace seqan3::detail

namespace std
{
//!\cond
template <typename tuple_t>
//!\cond
    requires (seqan3::detail::affine_score_cell<tuple_t> || seqan3::detail::affine_score_and_trace_cell<tuple_t>)
//!\endcond
struct tuple_size<seqan3::detail::affine_cell_proxy<tuple_t>> : public tuple_size<tuple_t>
{};

template <size_t index, typename tuple_t>
//!\cond
    requires (seqan3::detail::affine_score_cell<tuple_t> || seqan3::detail::affine_score_and_trace_cell<tuple_t>)
//!\endcond
struct tuple_element<index, seqan3::detail::affine_cell_proxy<tuple_t>> : public tuple_element<index, tuple_t>
{};
//!\endcond
} // namespace std
