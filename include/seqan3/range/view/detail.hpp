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
#include <range/v3/view/transform.hpp>

#include <seqan3/core/metafunction/range.hpp>
#include <seqan3/core/metafunction/template_inspection.hpp>
#include <seqan3/std/concept/range.hpp>

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
     * ```cpp
     * template <typename urng_t>
     *     requires seqan3::input_range_concept<urng_t> // or more/stricter requirements
     * struct view_foo
     * {
     *     // your implementation
     * };
     *
     * using foo_fn = generic_pipable_view_adaptor<view_foo>;
     *
     * namespace view
     * {
     * inline constexpr foo_fn foo;
     * }
     * ```
     *
     * This operator enables the function call syntax of the adaptor:
     *
     * ```cpp
     * std::vector<int> container{1, 2, 3};
     *
     * auto v = view::foo(container);    // if the view takes no constructor args beyond urange
     * auto v = view::foo(container, 7); // if the view takes e.g. an extra int argument
     * // in both cases v is now of type view_foo<std::vector<int>>
     * ```
     */
    template <typename urng_t, typename ...arg_types>
    auto operator()(urng_t && urange, arg_types && ...args) const
        requires input_range_concept<urng_t>
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
    struct auxiliary_functor_t
    {
        /*!\brief Define the operator() that is ultimately called inside the pipe-operator to fill in the `urange`
         * argument.
         * \tparam urng_t Type of the underlying range.
         * \param[in] urange The underlying range.
         */
        template <typename urng_t>
        auto operator()(urng_t && urange)
            requires (sizeof...(arg_types) > 0) && input_range_concept<urng_t>
        {
            return explode(std::forward<urng_t>(urange), std::make_index_sequence<sizeof...(arg_types)>{});
        }

        //!\brief Store the arguments.
        std::tuple<arg_types...> _arguments;

    protected:
        //!\brief Helper function to unpack the tuple.
        template <typename urng_t, size_t... Is>
        auto explode(urng_t && urange, std::index_sequence<Is...> const &)
        {
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
     * ```cpp
     * template <typename urng_t>
     *     requires seqan3::input_range_concept<urng_t> // or more/stricter requirements
     * struct view_foo
     * {
     *     // your implementation
     * };
     *
     * using foo_fn = generic_pipable_view_adaptor<view_foo>;
     *
     * namespace view
     * {
     * inline constexpr foo_fn foo;
     * }
     * ```
     *
     * This functor enables the function call syntax, without the actual underlying range:
     *
     * ```cpp
     * std::vector<int> container{1, 2, 3};
     *
     * auto v = view::foo(7); // v is NOT OF TYPE view_foo<std::vector<int>>
     *
     * // it's usually not used like above, instead use it inside a pipe:
     * auto v = container | view::foo(7);
     * ```
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
     * ```cpp
     * template <typename urng_t>
     *     requires seqan3::input_range_concept<urng_t> // or more/stricter requirements
     * struct view_foo
     * {
     *     // your implementation
     * };
     *
     * using foo_fn = generic_pipable_view_adaptor<view_foo>;
     *
     * namespace view
     * {
     * inline constexpr foo_fn foo;
     * }
     * ```
     *
     * This functor enables the pipe syntax:
     *
     * ```cpp
     * std::vector<int> container{1, 2, 3};
     *
     * auto v = container | view::foo(7);
     * //                 ^           ^
     * //     this operator           the intermediate operator() that returns a bound functor
     * ```
     */
    template <typename urng_t,
              typename ...arg_types>
    friend auto operator|(urng_t && urange,
                          auxiliary_functor_t<arg_types...> && bound_functor)
        requires (sizeof...(arg_types) > 0) && input_range_concept<urng_t>
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
     * ```cpp
     * template <typename urng_t>
     *     requires seqan3::input_range_concept<urng_t> // or more/stricter requirements
     * struct view_foo
     * {
     *     // your implementation
     * };
     *
     * using foo_fn = generic_pipable_view_adaptor<view_foo>;
     *
     * namespace view
     * {
     * inline constexpr foo_fn foo;
     * }
     * ```
     *
     * This functor enables the pipe syntax:
     *
     * ```cpp
     * std::vector<int> container{1, 2, 3};
     *
     * auto v = container | view::foo; // v is now of type view_foo<std::vector<int>
     * ```
     */
    template <typename urng_t>
    friend auto operator|(urng_t && urange,
                          derived_type const & fn)
        requires input_range_concept<urng_t>
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
 * ```cpp
 * // 1.
 * template <typename urng_t>
 *     requires seqan3::input_range_concept<urng_t> // or more/stricter requirements
 * struct view_foo
 * {
 *     // your implementation
 * };
 *
 * // 2. this is sufficient to declare the adaptor type:
 * using foo_fn = detail::generic_pipable_view_adaptor<view_foo>;
 *
 * namespace view
 * {
 * // 3. the adaptor instance
 * inline constexpr foo_fn foo;
 * }
 * ```
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
