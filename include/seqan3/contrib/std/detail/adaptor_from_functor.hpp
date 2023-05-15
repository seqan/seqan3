// -----------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2023, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2023, Knut Reinert & MPI für molekulare Genetik
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at: https://github.com/seqan/seqan-std/blob/main/LICENSE
// -----------------------------------------------------------------------------------------------------

/*!\file
 * \brief Provides seqan::std::detail::adaptor_from_functor
 * \author Hannes Hauswedell <hannes.hauswedell AT fu-berlin.de>
 */

// File might be included from multiple libraries.
#ifndef SEQAN_STD_DETAIL_ADAPTOR_FROM_FUNCTOR
#define SEQAN_STD_DETAIL_ADAPTOR_FROM_FUNCTOR

#include "adaptor_base.hpp"

namespace seqan::std::detail
{

// ============================================================================
//  adaptor_from_functor
// ============================================================================

//!\brief Template for range adaptor closure objects that store arguments and wrap a proto-adaptor.
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
    template <::std::ranges::input_range urng_t>
    constexpr auto impl(urng_t && urange, stored_args_ts... args) const
    {
        return fun(::std::forward<urng_t>(urange), ::std::forward<stored_args_ts>(args)...);
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
        base_type{::std::forward<stored_args_ts>(args)...},
        fun{::std::move(f)}
    {}
    //!\}
};

} // namespace seqan::std::detail

#endif // SEQAN_STD_DETAIL_ADAPTOR_FROM_FUNCTOR
