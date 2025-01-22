// SPDX-FileCopyrightText: 2006-2025 Knut Reinert & Freie Universität Berlin
// SPDX-FileCopyrightText: 2016-2025 Knut Reinert & MPI für molekulare Genetik
// SPDX-License-Identifier: BSD-3-Clause

/*!\file
 * \brief Provides seqan3::detail::adaptor_from_functor
 * \author Hannes Hauswedell <hannes.hauswedell AT fu-berlin.de>
 */

#pragma once

#include <seqan3/core/range/detail/adaptor_base.hpp>

namespace seqan3::detail
{

// ============================================================================
//  adaptor_from_functor
// ============================================================================

/*!\brief Template for range adaptor closure objects that store arguments and wrap a proto-adaptor.
 * \tparam functor_type Type of the proto-adaptor functor.
 * \tparam stored_args_ts Types of the stored arguments.
 * \ingroup core_range
 *
 * \details
 *
 * This template is particularly useful in combination with range adaptor objects that are not closure
 * objects (a.k.a. proto adaptors).
 * These objects take other parameters in addition to the range parameter that they might need to store
 * in which case they can return an object of this type with the stored arguments.
 *
 * This template delegates back to the proto adaptor for handling range input (impl() calls operator() of the
 * proto-adaptor).
 *
 * # Example
 *
 * From include/seqan3/search/views/kmer_hash.hpp:
 *
 * \snippet include/seqan3/search/views/kmer_hash.hpp adaptor_def
 *
 * This is the full proto-adaptor, first look at the second member function: it handles range and argument input and
 * delegates to the view's constructor. In other, simpler cases you could invoke other adaptors here.
 *
 * And it provides an `operator()` that takes only the argument and returns a range adaptor closure object
 * with the argument wrapped (first member function). The proto-adaptor is passed to the closure object, as well, so
 * that the proto-adaptor's `operator()` can be used for handling range input.
 *
 * Note that the proto-adaptor does not provide `operator|`, this is only required of adaptor closure objects.
 */
template <typename functor_type, typename... stored_args_ts>
class adaptor_from_functor :
    public adaptor_base<adaptor_from_functor<functor_type, stored_args_ts...>, stored_args_ts...>
{
private:
    //!\brief Type of the CRTP-base.
    using base_type = adaptor_base<adaptor_from_functor<functor_type, stored_args_ts...>, stored_args_ts...>;

    //!\brief Befriend the base class so it can call impl().
    friend base_type;

    //!\brief The stored functor, usually a "proto-adaptor".
    functor_type fun;

    /*!\brief Delegate the range argument and stored arguments to the wrapped functor.
     * \tparam    urng_t         The underlying range type.
     * \tparam    stored_args_ts The arguments to the view (this first one will be a range, the rest is optional).
     * \param[in] urange         The underling range.
     * \param[in] args           The arguments to the constructor.
     * \returns Whatever the wrapped functor returns, usually a view.
     */
    template <std::ranges::input_range urng_t>
    constexpr auto impl(urng_t && urange, stored_args_ts... args) const
    {
        return fun(std::forward<urng_t>(urange), std::forward<stored_args_ts>(args)...);
    }

public:
    /*!\name Constructors, destructor and assignment
     * \{
     */
    constexpr adaptor_from_functor() = default;                                                  //!< Defaulted.
    constexpr adaptor_from_functor(adaptor_from_functor const &) noexcept = default;             //!< Defaulted.
    constexpr adaptor_from_functor(adaptor_from_functor &&) noexcept = default;                  //!< Defaulted.
    constexpr adaptor_from_functor & operator=(adaptor_from_functor const &) noexcept = default; //!< Defaulted.
    constexpr adaptor_from_functor & operator=(adaptor_from_functor &&) noexcept = default;      //!< Defaulted.
    ~adaptor_from_functor() noexcept = default;                                                  //!< Defaulted.

    //!\brief Construct from functor and possibly arguments.
    constexpr adaptor_from_functor(functor_type f, stored_args_ts... args) :
        base_type{std::forward<stored_args_ts>(args)...},
        fun{std::move(f)}
    {}
    //!\}
};

} // namespace seqan3::detail
