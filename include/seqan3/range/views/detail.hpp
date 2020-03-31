// -----------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2020, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2020, Knut Reinert & MPI für molekulare Genetik
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at: https://github.com/seqan/seqan3/blob/master/LICENSE.md
// -----------------------------------------------------------------------------------------------------

/*!\file
 * \ingroup views
 * \brief Auxiliary header for the \link views views submodule \endlink.
 * \author Hannes Hauswedell <hannes.hauswedell AT fu-berlin.de>
 */

#pragma once

#include <tuple>

#include <seqan3/core/type_traits/basic.hpp>
#include <seqan3/core/type_traits/range.hpp>
#include <seqan3/std/ranges>

namespace seqan3::detail
{

// ============================================================================
//  forwards
// ============================================================================

template <typename left_adaptor_t, typename right_adaptor_t>
class combined_adaptor;

// ============================================================================
//  adaptor_base
// ============================================================================

/*!\brief CRTP-base to simplify the definition of range adaptor closure objects and similar types.
 * \tparam derived_type CRTP specialisation parameter.
 * \tparam stored_args_ts Some closure objects store values, these are their types.
 * \ingroup views
 *
 * \details
 *
 * See \ref howto_write_a_view for details on what range adaptor closure objects are.
 *
 * This class requires from the `derived_type`:
 *   * an `impl()` function that takes a range and possibly further arguments.
 *
 * This base class provides:
 *   * Storage of args, including perfect forwarding and moving-out.
 *   * `operator()` that takes a range and passes the range and stored args to the derived_type.
 *   * `operator|` that also takes a range and passes the range and stored args to the derived_type.
 *   * `operator|` that takes another adaptor and returns seqan3::detail::combined_adaptor
 *
 * # Examples
 *
 * ## "Direct derivation"
 *
 * From include/seqan3/range/views/persist.hpp:
 *
 * \snippet include/seqan3/range/views/persist.hpp adaptor_def
 *
 * This adaptor directly derives from adaptor_base (instead of just using
 * `seqan3::detail::adaptor_for_view_without_args`) so that it can decide between
 * actually calling the respective view's constructor and delegating to a different view
 * (in this case `std::views::all`) based on the type of input.
 *
 * ## Derived templates
 *
 * * `seqan3::detail::combined_adaptor` : Combine two adaptor closure objects into a new one.
 * * `seqan3::detail::adaptor_for_view_without_args` : Create an adaptor closure object for a view that doesn't require
 *   args.
 * * `seqan3::detail::adaptor_from_functor` : Create an adaptor closure object with possibly stored arguments that
 *   delegates to a functor, usually a stored proto-adaptor (for view classes that require arguments or adaptors that
 *   just forward to other adaptors).
 * * `seqan3::views::deep` : Wraps an adaptor closure object or proto adaptor object and modifies the behaviour.
 */
template <typename derived_type, typename ...stored_args_ts>
class adaptor_base
{
private:
    //!\brief Stores the arguments.
    std::tuple<stored_args_ts...> arguments;

    //!\brief Helper function to unpack the tuple and delegate to the derived type.
    template <typename urng_t, size_t... Is>
    constexpr auto pass_args_to_impl(urng_t && urange, std::index_sequence<Is...> const &) const &
    {
        // std::get returns lvalue-reference to value, but we need to copy the values
        // so that the view does not depend on the functor
        return static_cast<derived_type const &>(*this).impl(
            std::forward<urng_t>(urange),
            std::tuple_element_t<Is, std::tuple<stored_args_ts...>>(std::get<Is>(arguments))...);
    }

    //!\overload
    template <typename urng_t, size_t... Is>
    constexpr auto pass_args_to_impl(urng_t && urange, std::index_sequence<Is...> const &) &&
    {
        // move out values, because we don't need them anymore (*this is temporary)
        return static_cast<derived_type &&>(*this).impl(
            std::forward<urng_t>(urange),
            std::tuple_element_t<Is, std::tuple<stored_args_ts...>>(std::get<Is>(std::move(arguments)))...);
    }


    //!\brief Befriend the derived_type so it can access private members if need be.
    friend derived_type;

public:
    /*!\name Constructors, destructor and assignment
     * \{
     */
    // a default constructor is not provided, however the constructor below might be one.
    constexpr adaptor_base(adaptor_base const &)                noexcept = default; //!< Defaulted.
    constexpr adaptor_base(adaptor_base &&)                     noexcept = default; //!< Defaulted.
    constexpr adaptor_base & operator=(adaptor_base const &)    noexcept = default; //!< Defaulted.
    constexpr adaptor_base & operator=(adaptor_base &&)         noexcept = default; //!< Defaulted.
    ~adaptor_base()                                             noexcept = default; //!< Defaulted.

    //!\brief Constructor with possible arguments; becomes a default constructor for adaptors without args.
    constexpr adaptor_base(stored_args_ts ... args)
        noexcept(noexcept(std::tuple<stored_args_ts...>{std::forward<stored_args_ts>(args)...})) :
        arguments{std::forward<stored_args_ts>(args)...}
    {}
    //!\}

    //!\brief Function-style overload for ranges.
    template <std::ranges::input_range urng_t>
    constexpr auto operator()(urng_t && urange) const &
    {
        return pass_args_to_impl(std::forward<urng_t>(urange), std::make_index_sequence<sizeof...(stored_args_ts)>{});
    }

    //!\overload
    template <std::ranges::input_range urng_t>
    constexpr auto operator()(urng_t && urange) &&
    {
        return std::move(*this).pass_args_to_impl(std::forward<urng_t>(urange),
                                                  std::make_index_sequence<sizeof...(stored_args_ts)>{});
    }

    /*!\brief Left-hand-side pipe operator that handles range input or composing.
     * \tparam    arg_t LHS argument type.
     * \param[in] arg   LHS, either a range or another adaptor.
     * \param[in] me    RHS.
     *
     * \details
     *
     * If the LHS models std::ranges::input_range this functions delegates to operator()
     * otherwise it assumes the LHS is another range adaptor closure object and it
     * returns a seqan3::detail::combined_adaptor that wraps that adaptor and this one.
     */
    template <typename arg_t>
    constexpr friend auto operator|(arg_t && arg, derived_type const & me)
    {
        if constexpr (std::ranges::input_range<arg_t>)
            return me(std::forward<arg_t>(arg));
        else
            return combined_adaptor{std::forward<arg_t>(arg), me};
    }

    //!\overload
    template <typename arg_t>
    constexpr friend auto operator|(arg_t && arg, derived_type && me)
    {
        if constexpr (std::ranges::input_range<arg_t>)
            return std::move(me)(std::forward<arg_t>(arg));
        else
            return combined_adaptor{std::forward<arg_t>(arg), std::move(me)};
    }

    /*!\brief Right-hand-side pipe operator that handles composing.
     * \tparam arg_t RHS argument type.
     * \param[in] me  LHS, another adaptor.
     * \param[in] arg RHS.
     *
     * \details
     *
     * This operator assumes that the RHS is another range adaptor closure object and it
     * returns a seqan3::detail::combined_adaptor that wraps that adaptor and this one.
     *
     * Note that this operator takes adaptor_base and not derived_type as parameter type
     * so that an implicit conversion is required for invoking this function.
     * This prioritises the LHS-pipe-operator overload and resolves ambiguity in cases
     * where they are competing.
     * It has no semantic meaning, the type is down-cast inside the function again.
     */
    template <typename arg_t>
    constexpr friend auto operator|(adaptor_base const & me, arg_t && arg)
    {
        return combined_adaptor{static_cast<derived_type const &>(me), std::forward<arg_t>(arg)};
    }

    //!\overload
    template <typename arg_t>
    constexpr friend auto operator|(adaptor_base && me, arg_t && arg)
    {
        return combined_adaptor{static_cast<derived_type &&>(me), std::forward<arg_t>(arg)};
    }
};

// ============================================================================
//  combined_adaptor
// ============================================================================

/*!\brief Template for range adaptor closure objects that consist of two other range adaptor closure objects.
 * \tparam left_adaptor_t  Type of the first stored adaptor.
 * \tparam right_adaptor_t Type of the second stored adaptor.
 * \ingroup views
 *
 * \details
 *
 * If invoked with a range, this adaptor resolves to piping the range into left adaptor and the resulting
 * range into the right adaptor.
 */
template <typename left_adaptor_t, typename right_adaptor_t>
class combined_adaptor :
    public adaptor_base<combined_adaptor<left_adaptor_t, right_adaptor_t>,
                        left_adaptor_t,
                        right_adaptor_t>
{
private:
    //!\brief Type of the CRTP-base.
    using base_type = adaptor_base<combined_adaptor<left_adaptor_t, right_adaptor_t>,
                                   left_adaptor_t,
                                   right_adaptor_t>;

    //!\brief Befriend the base class so it can call impl().
    friend base_type;

    //!\brief Combine all arguments via `operator|`.
    template <std::ranges::input_range urng_t,
              typename left_adaptor_t_,
              typename right_adaptor_t_>
    static auto impl(urng_t && urange, left_adaptor_t_ && left_adaptor, right_adaptor_t_ && right_adaptor)
    {
        return std::forward<urng_t>(urange)
             | std::forward<left_adaptor_t_>(left_adaptor)
             | std::forward<right_adaptor_t_>(right_adaptor);
    }

public:
    /*!\name Constructors, destructor and assignment
     * \{
     */
    constexpr combined_adaptor()                                              = default; //!< Defaulted.
    constexpr combined_adaptor(combined_adaptor const &)             noexcept = default; //!< Defaulted.
    constexpr combined_adaptor(combined_adaptor &&)                  noexcept = default; //!< Defaulted.
    constexpr combined_adaptor & operator=(combined_adaptor const &) noexcept = default; //!< Defaulted.
    constexpr combined_adaptor & operator=(combined_adaptor &&)      noexcept = default; //!< Defaulted.
    ~combined_adaptor()                                              noexcept = default; //!< Defaulted.

    //!\brief Inherit the base type's constructors.
    using base_type::base_type;

    //!\brief Store both arguments in the adaptor.
    constexpr combined_adaptor(left_adaptor_t l, right_adaptor_t r) :
        base_type{std::forward<left_adaptor_t>(l), std::forward<right_adaptor_t>(r)}
    {}
    //!\}
};

// ============================================================================
//  adaptor_for_view_without_args
// ============================================================================

/*!\brief Template for range adaptor closure objects that store no arguments and always delegate to the
 *        view constructor.
 * \tparam view_type The view type template.
 * \ingroup views
 *
 * \details
 *
 * Use this adaptor template when you always want to delegate to the view's constructor and you have no arguments.
 * Since it's a one-line it's easier than specialising seqan3::detail::adaptor_base.
 *
 * # Example
 *
 * (from include/seqan3/range/views/single_pass_input.hpp)
 *
 * This is the signature of the view type template in namespace seqan3::detail:
 *
 * \snippet include/seqan3/range/views/single_pass_input.hpp view_def
 *
 * This is the definition of the range adaptor closure object, it will always delegate to the
 * constructor: seqan3::detail::single_pass_input_view::single_pass_input_view():
 *
 * \snippet include/seqan3/range/views/single_pass_input.hpp adaptor_def
 */
template <template <typename, typename ...> typename view_type>
class adaptor_for_view_without_args : public adaptor_base<adaptor_for_view_without_args<view_type>>
{
private:
    //!\brief Type of the CRTP-base.
    using base_type = adaptor_base<adaptor_for_view_without_args<view_type>>;

    //!\brief Befriend the base class so it can call impl().
    friend base_type;

    /*!\brief Call the view's constructor with the given arguments (all of the base class'es operators ultimately
     * resolve to this function call).
     * \tparam    arg_types The arguments to the view (this first one will be a range, the rest is optional).
     * \param[in] args      The arguments to the constructor.
     * \returns An instance of `view_type`.
     */
    template <typename ...arg_types>
    static auto impl(arg_types && ... args)
    {
        return view_type{std::forward<arg_types>(args)...};
    }

public:
    /*!\name Constructors, destructor and assignment
     * \{
     */
    //!\brief Defaulted.
    constexpr adaptor_for_view_without_args()                                                           = default;
    //!\brief Defaulted.
    constexpr adaptor_for_view_without_args(adaptor_for_view_without_args const &)             noexcept = default;
    //!\brief Defaulted.
    constexpr adaptor_for_view_without_args(adaptor_for_view_without_args &&)                  noexcept = default;
    //!\brief Defaulted.
    constexpr adaptor_for_view_without_args & operator=(adaptor_for_view_without_args const &) noexcept = default;
    //!\brief Defaulted.
    constexpr adaptor_for_view_without_args & operator=(adaptor_for_view_without_args &&)      noexcept = default;
    //!\brief Defaulted.
    ~adaptor_for_view_without_args()                                                           noexcept = default;

    //!\brief Inherit the base type's constructors.
    using base_type::base_type;
    //!\}
};

// ============================================================================
//  adaptor_from_functor
// ============================================================================

/*!\brief Template for range adaptor closure objects that store arguments and wrap a proto-adaptor.
 * \tparam functor_type Type of the proto-adaptor functor.
 * \tparam stored_args_ts Types of the stored arguments.
 * \ingroup views
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
 * From include/seqan3/range/views/kmer_hash.hpp:
 *
 * \snippet include/seqan3/range/views/kmer_hash.hpp adaptor_def
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
template <typename functor_type, typename ...stored_args_ts>
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
    constexpr auto impl(urng_t && urange, stored_args_ts ... args) const
    {
        return fun(std::forward<urng_t>(urange), std::forward<stored_args_ts>(args)...);
    }

public:
    /*!\name Constructors, destructor and assignment
     * \{
     */
    constexpr adaptor_from_functor()                                                  = default; //!< Defaulted.
    constexpr adaptor_from_functor(adaptor_from_functor const &)             noexcept = default; //!< Defaulted.
    constexpr adaptor_from_functor(adaptor_from_functor &&)                  noexcept = default; //!< Defaulted.
    constexpr adaptor_from_functor & operator=(adaptor_from_functor const &) noexcept = default; //!< Defaulted.
    constexpr adaptor_from_functor & operator=(adaptor_from_functor &&)      noexcept = default; //!< Defaulted.
    ~adaptor_from_functor()                                                  noexcept = default; //!< Defaulted.

    //!\brief Construct from functor and possibly arguments.
    constexpr adaptor_from_functor(functor_type f, stored_args_ts ... args) :
        base_type{std::forward<stored_args_ts>(args)...}, fun{std::move(f)}
    {}
    //!\}
};

} // namespace seqan3::detail
