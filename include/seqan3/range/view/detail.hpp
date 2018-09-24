// ============================================================================
//                 SeqAn - The Library for Sequence Analysis
// ============================================================================
//
// Copyright (c) 2006-2018, Knut Reinert & Freie Universitaet Berlin
// Copyright (c) 2016-2018, Knut Reinert & MPI Molekulare Genetik
// All rights reserved.
//
// Redistribution and use in source and binary forms, with or without
// modification, are permitted provided that the following conditions are met:
//
//     * Redistributions of source code must retain the above copyright
//       notice, this list of conditions and the following disclaimer.
//     * Redistributions in binary form must reproduce the above copyright
//       notice, this list of conditions and the following disclaimer in the
//       documentation and/or other materials provided with the distribution.
//     * Neither the name of Knut Reinert or the FU Berlin nor the names of
//       its contributors may be used to endorse or promote products derived
//       from this software without specific prior written permission.
//
// THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
// AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
// IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
// ARE DISCLAIMED. IN NO EVENT SHALL KNUT REINERT OR THE FU BERLIN BE LIABLE
// FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL
// DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR
// SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER
// CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT
// LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY
// OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH
// DAMAGE.
//
// ============================================================================

/*!\file
 * \ingroup view
 * \brief Auxiliary header for the \link view view submodule \endlink.
 * \author Hannes Hauswedell <hannes.hauswedell AT fu-berlin.de>
 */

#pragma once

#include <range/v3/range_fwd.hpp>

#include <seqan3/core/metafunction/range.hpp>
#include <seqan3/core/metafunction/template_inspection.hpp>
#include <seqan3/core/type_list.hpp>
#include <seqan3/std/ranges>

namespace seqan3::detail
{

// ============================================================================
//  view_base
// ============================================================================

/*!\brief An empty base class that our views inherit from so they are treated as views by the STL/range-v3.
 * \ingroup view
 */
struct view_base : public ranges::view_base
{};

// ============================================================================
//  size_type_t_or_void
// ============================================================================

//!\brief Transformation trait that resolves to seqan3::size_type_t if possible and void otherwise. [Default Overload]
template <typename t>
struct size_type_t_or_void
{
    //!\brief The default is to expose `void`.
    using type = void;
};

//!\brief Transformation trait that resolves to seqan3::size_type_t if possible and void otherwise. [Specialisation]
template <std::ranges::SizedRange t>
struct size_type_t_or_void<t>
{
    //!\brief For sized ranges the seqan3::size_type_t is exposed.
    using type = size_type_t<t>;
};

// ============================================================================
//  pipable_adaptor_base
// ============================================================================

/*!\brief   A CRTP base class for implementing the pipable behaviour of views.
 * \ingroup view
 * \tparam  derived_type The CRTP specialisation, must provide a `impl()`-member.
 *
 * \details
 *
 * ### Background
 *
 * A full view implementation consists of three entities:
 *
 *   1. the actual view, e.g. `view_foo`;
 *   2. an adaptor that returns instances of 1. and enables usage in pipes, e.g. `foo_fn`
 *   3. an instance of the adaptor, e.g. `view::foo` that is usable by consumers of your view
 *
 * For more details, see the HowTo on writing views in the SeqAn3 wiki: https://github.com/seqan/seqan3/wiki
 *
 * ### Details
 *
 * This CRTP base class can be used to simplify the implementation of the adaptor type. For most use cases it is
 * sufficient to just use the specialisation seqan3::detail::generic_pipable_view_adaptor, but other specialisations
 * like seqan3::detail::deep_view_adaptor_wrapper also make use of this base template.
 */
template <typename derived_type>
class pipable_adaptor_base
{
public:
    /*!\brief       Enables function style usage of the adaptor, forwards to the actual implementation in `derived_type`.
     * \tparam      urng_t Type of the underlying range.
     * \tparam      arg_types The types of further arguments.
     * \param[in]   urange The underlying range.
     * \param[in]   args Further arguments.
     * \returns     Depends on the `derived_type`.
     *
     * \details
     *
     * This function works for views that take arguments, and for view that take no arguments, i.e.
     * depending on the view args needs to be exactly those parameters that view expects which
     * can/must be none for certain views.
     *
     * ### Example (based on the derived type seqan3::detail::generic_pipable_view_adaptor)
     *
     * Given:
     *
     * \snippet test/snippet/range/view/detail.cpp usage
     *
     * This operator enables the function call syntax of the adaptor:
     *
     * \snippet test/snippet/range/view/detail.cpp function_call
     */
    template <std::ranges::InputRange urng_t, typename ...arg_types>
    auto operator()(urng_t && urange, arg_types && ...args) const
    {
        return static_cast<derived_type const *>(this)->impl(std::forward<urng_t>(urange),
                                                             std::forward<arg_types>(args)...);
    }

    /*!\brief   Binds operator() with a placeholder for range so that the operator() without `urange` argument works
     *          [for adaptors with arguments].
     * \tparam  arg_types The types of the arguments.
     *
     * \details
     *
     * This replaces std::bind which doesn't work here because of the parent being template template and thus
     * the type not being fully resolved. A lambda function wrapper could be used, but its type wouldn't be
     * deducible which we need for the definition of the pipe operator.
     *
     * TODO re-evaluate this since pipable_adaptor_base is no longer template template
     */
    template <typename ...arg_types>
    //!\cond
        requires (sizeof...(arg_types) > 0)
    //!\endcond
    struct auxiliary_functor_t
    {
        /*!\brief Define the operator() that is ultimately called inside the pipe-operator to fill in the `urange`
         * argument.
         * \tparam urng_t Type of the underlying range.
         * \param[in] urange The underlying range.
         */
        template <std::ranges::InputRange urng_t>
        auto operator()(urng_t && urange)
        {
            return explode(std::forward<urng_t>(urange), std::make_index_sequence<sizeof...(arg_types)>{});
        }

        //!\overload
        template <std::ranges::InputRange urng_t>
        auto operator()(urng_t && urange) const
        {
            return explode(std::forward<urng_t>(urange), std::make_index_sequence<sizeof...(arg_types)>{});
        }

        //!\brief Store the arguments.
        std::tuple<arg_types...> _arguments;

    protected:
        //!\brief Helper function to unpack the tuple.
        template <typename urng_t, size_t... Is>
        auto explode(urng_t && urange, std::index_sequence<Is...> const &) const &
        {
            // std::get returns lvalue-reference to value, but we need to copy the values
            return pipable_adaptor_base{}(
                std::forward<urng_t>(urange),
                static_cast<std::tuple_element_t<Is, std::tuple<arg_types...>>>(std::get<Is>(_arguments))...);
        }

        //!\brief Helper function to unpack the tuple.
        template <typename urng_t, size_t... Is>
        auto explode(urng_t && urange, std::index_sequence<Is...> const &) &
        {
            // std::get returns lvalue-reference to value, but we need to copy the values
            return pipable_adaptor_base{}(
                std::forward<urng_t>(urange),
                static_cast<std::tuple_element_t<Is, std::tuple<arg_types...>>>(std::get<Is>(_arguments))...);
        }

        //!\overload
        template <typename urng_t, size_t... Is>
        auto explode(urng_t && urange, std::index_sequence<Is...> const &) &&
        {
            // move out values, because we don't need them anymore (this is temporary)
            return pipable_adaptor_base{}(std::forward<urng_t>(urange),
                                          std::forward<arg_types>(std::get<Is>(_arguments))...);
        }
    };

    /*!\brief       Enables calling the functor without a `urange` parameter, necessary inside a pipe statement
     *              [for adaptors with arguments].
     * \tparam      arg_types The types of the arguments.
     * \param[in]   args Arguments stored in the auxiliary functor.
     * \returns     An instance of auxiliary_functor_t for use in the pipe notation, **not an instance of `view_type`**!
     *
     * \details
     *
     * ### Example (based on the derived type seqan3::detail::generic_pipable_view_adaptor)
     *
     * Given:
     *
     * \snippet test/snippet/range/view/detail.cpp usage
     *
     * This functor enables the function call syntax, without the actual underlying range:
     *
     * \snippet test/snippet/range/view/detail.cpp function_call_2
     */
    template <typename ...arg_types>
    auxiliary_functor_t<arg_types...> operator()(arg_types && ...args) const
    {
        return {{std::forward<arg_types>(args)...}};
    }

    /*!\brief       Enables pipe style usage of the view, forwards indirectly to the constructor of `view_type`
     *              [for adaptors with arguments].
     * \tparam      urng_t Type of the underlying range.
     * \tparam      arg_types The types of the arguments.
     * \param[in]   urange The underlying range [left hand side argument to `|`].
     * \param[in]   bound_functor A functor of type auxiliary_functor_t, the result of a call to operator()()
     *              [right hand side argument to `|`].
     * \returns     An object of type `view_type`.
     *
     * \details
     *
     * ### Example (based on the derived type seqan3::detail::generic_pipable_view_adaptor)
     *
     * Given:
     *
     * \snippet test/snippet/range/view/detail.cpp usage
     *
     * This functor enables the pipe syntax:
     *
     * \snippet test/snippet/range/view/detail.cpp pipe_syntax
     */
    template <std::ranges::InputRange urng_t, typename ...arg_types>
    friend auto operator|(urng_t && urange, auxiliary_functor_t<arg_types...> & bound_functor)
        requires (sizeof...(arg_types) > 0)
    {
        return bound_functor(std::forward<urng_t>(urange));
    }

    //!\overload
    template <std::ranges::InputRange urng_t, typename ...arg_types>
    friend auto operator|(urng_t && urange, auxiliary_functor_t<arg_types...> const & bound_functor)
        requires (sizeof...(arg_types) > 0)
    {
        return bound_functor(std::forward<urng_t>(urange));
    }

    //!\overload
    template <std::ranges::InputRange urng_t, typename ...arg_types>
    friend auto operator|(urng_t && urange, auxiliary_functor_t<arg_types...> && bound_functor)
        requires (sizeof...(arg_types) > 0)
    {
        return bound_functor(std::forward<urng_t>(urange));
    }

    /*!\brief       Enables pipe style usage of the view, forwards to the constructor of view_type
     *              [for adaptors without arguments].
     * \tparam      urng_t Type of the underlying range.
     * \param[in]   urange The underlying range [left hand side argument to `|`].
     * \param[in]   fn An instance of the functor, e.g. `view::foo` [right hand side argument to `|`].
     * \returns     An object of type `view_type`.
     *
     * \details
     *
     * ### Example (based on the derived type seqan3::detail::generic_pipable_view_adaptor)
     *
     * Given:
     *
     * \snippet test/snippet/range/view/detail.cpp usage
     *
     * This functor enables the pipe syntax:
     *
     * \snippet test/snippet/range/view/detail.cpp pipe_syntax_2
     */
    template <std::ranges::InputRange urng_t>
    friend auto operator|(urng_t && urange,
                          derived_type const & fn)
    {
        return fn(std::forward<urng_t>(urange));
    }

private:
    /*!\name Constructors, destructor and assignment
     * \brief All are declared private to prevent direct use of the CRTP base.
     * \sa https://isocpp.org/blog/2017/04/quick-q-prevent-user-from-derive-from-incorrect-crtp-base
     * \{
     */
    pipable_adaptor_base() = default;
    constexpr pipable_adaptor_base(pipable_adaptor_base const &) = default;
    constexpr pipable_adaptor_base(pipable_adaptor_base &&) = default;
    constexpr pipable_adaptor_base & operator =(pipable_adaptor_base const &) = default;
    constexpr pipable_adaptor_base & operator =(pipable_adaptor_base &&) = default;
    ~pipable_adaptor_base() = default;
    //!\brief befriend the derived type so that it can instantiate
    friend derived_type;
    //!\}
};

// ============================================================================
//  generic_pipable_view_adaptor
// ============================================================================

/*!\brief   A class template template that makes a view pipeable.
 * \ingroup view
 * \tparam  view_type The view template that your wish to define the adaptor for.
 *
 * \details
 *
 * A full view implementation consists of three entities:
 *
 *   1. the actual view, e.g. `view_foo`;
 *   2. an adaptor that returns instances of 1. and enables usage in pipes, e.g. `foo_fn`
 *   3. an instance of the adaptor, e.g. `view::foo` that is usable by consumers of your view
 *
 * For more details, see the HowTo on writing views in the SeqAn3 wiki: https://github.com/seqan/seqan3/wiki
 *
 * This template generates the adaptor for you:
 *
 * \snippet test/snippet/range/view/detail.cpp adaptor_template
 */
template <template <typename, typename...> class view_type>
class generic_pipable_view_adaptor : public pipable_adaptor_base<generic_pipable_view_adaptor<view_type>>
{
private:
    //!\brief Type of the CRTP-base.
    using base_type = pipable_adaptor_base<generic_pipable_view_adaptor<view_type>>;

public:
    //!\brief Inherit the base class's Constructors.
    using base_type::base_type;

private:
    //!\brief Befriend the base class so it can call impl().
    friend base_type;

    /*!\brief       Call the view's constructor with the given arguments (all of the base class'es operators ultimately
     *              resolve to this function call).
     * \tparam      arg_types The arguments to the view (this first one will be a range, the rest is optional).
     * \param[in]   args The arguments to the constructor.
     * \returns     An instance of `view_type`.
     */
    template <typename ... arg_types>
    static auto impl(arg_types && ... args)
    {
        return view_type{std::forward<arg_types>(args)...};
    }
};

} // namespace seqan3::detail
