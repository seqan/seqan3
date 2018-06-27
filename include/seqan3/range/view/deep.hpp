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
 * \author Hannes Hauswedell <hannes.hauswedell AT fu-berlin.de>
 * \brief Provides seqan3::view::deep.
 */

#pragma once

#include <range/v3/view/transform.hpp>

#include <seqan3/range/view/detail.hpp>

namespace seqan3::view
{

/*!\brief   A wrapper type around an existing view adaptor that enables "deep view" behaviour for that view.
 * \ingroup view
 * \tparam  underlying_adaptor_t The type of the adaptor being wrapped.
 *
 * \details
 *
 * ### Deep views
 *
 * If you pass a range to a view that view performs some transformation on that range. If the range passed is
 * multi-dimensional (i.e. a range-of-ranges) that transformation happens on the outermost range. So if you
 * call ranges::view::reverse on a range-of-dna-ranges, it will revert *the order* of the dna-ranges, but leave
 * the dna-ranges themselves unchanged.
 *
 * In some cases this is not desirable or even possible, i.e. seqan3::view::complement performs it's operation on
 * nucleotide-ranges and it would be logical to do so, even it is passed a range-of-nucleotide-ranges (it obviously
 * cannot transform the outer range). We call these views "deep views" as they always perform their operation on
 * the innermost ranges of a multi-dimensional range; in case the input is a one-dimensional range, deepness
 * does not modify the behaviour.
 *
 * ### Using view::deep
 *
 * Strictly speaking, seqan3::view::deep is a view adaptor adaptor, i.e. it gets passed **another adaptor when being
 * constructed** (not via the pipe!) and returns an adaptor that behaves like the underlying one, except being deep.
 *
 * You can use it mostly like any other view (adaptor) with some subtle differences, illustrated in the examples below.
 *
 * ### View properties
 *
 * The returned view has the same requirements and guarantees as those of the underlying adaptor type, except that
 * it is also deep, i.e. if the underlying range is range-of-ranges, all transformations apply to the innermost ranges
 * and conversely the requirements also apply to the innermost ranges of the underlying range and guarantees apply
 * to the innermost ranges of the returned range.
 *
 * *For the higher dimensions* (all except the innermost ranges) the following properties hold:
 *
 * | range concepts and reference_t      | `urng_t` (underlying range type)      | `rrng_t` (returned range type)                     |
 * |-------------------------------------|:-------------------------------------:|:--------------------------------------------------:|
 * | seqan3::input_range_concept         | *required*                            | *preserved*                                        |
 * | seqan3::forward_range_concept       |                                       | *preserved*                                        |
 * | seqan3::bidirectional_range_concept |                                       | *preserved*                                        |
 * | seqan3::random_access_range_concept |                                       | *preserved*                                        |
 * |                                     |                                       |                                                    |
 * | seqan3::view_concept                |                                       | *guaranteed*                                       |
 * | seqan3::sized_range_concept         |                                       | *preserved*                                        |
 * | seqan3::bounded_range_concept       |                                       | *preserved*                                        |
 * | seqan3::output_range_concept        |                                       | *lost*                                             |
 * | seqan3::const_iterable_concept      |                                       | *preserved*                                        |
 * |                                     |                                       |                                                    |
 * | seqan3::reference_t                 | seqan3::input_range_concept           | seqan3::input_range_concept + seqan3::view_concept |
 *
 * ### Examples
 *
 * Wrapping an adaptor that takes no parameters:
 *
 * ```cpp
 * std::vector<dna5_vector> foo{"AAATTT"_dna5, "CCCGGG"_dna5};
 *
 * auto r = foo | ranges::view::reverse;             // == [ [C,C,C,G,G,G], [A,A,A,T,T,T] ]
 *
 * auto d = foo | view::deep{ranges::view::reverse}; // == [ [T,T,T,A,A,A], [G,G,G,C,C,C] ]
 *
 * // You can also create a permanent alias:
 * namespace view
 * {
 * inline auto const deep_reverse = deep{ranges::view::reverse};
 * }
 *
 * auto e = foo | view::deep_reverse;                // == [ [T,T,T,A,A,A], [G,G,G,C,C,C] ]
 * ```
 *
 * Wrapping an adaptor that takes parameters:
 *
 * ```cpp
 *
 * std::vector<dna5_vector> foo{"AAATTT"_dna5, "CCCGGG"_dna5};
 *
 * auto t = foo | ranges::view::take(1);             // == [ [A,A,A,T,T,T] ]
 *
 * auto d = foo | view::deep{ranges::view::take}(1); // == [ [A], [C] ]
 * // constructor arguments passed via {} and arguments to underlying view passed via ()
 *
 * // In this case especially, an alias improves readability:
 * namespace view
 * {
 * inline auto const deep_take = deep{ranges::view::take};
 * }
 *
 * auto e = foo | view::deep_take(1);                // == [ [A], [C] ]
 * ```
 *
 * The above example illustrates that view::deep has two sets of arguments, the **arguments to construct** this adaptor,
 * and the arguments passed to the underlying adaptor when **calling** this adaptor. You can use `()` for both, but
 * we highly recommend to use `{}` to not confuse these; or just use an alias.
 *
 * \attention Note that in the case of parameter handling the arguments to view::deep are **copied** to each invocation
 * of the underlying adaptor if they are temporaries. This is no problem for small objects like the integer above,
 * but might be expensive for larger ones. To avoid this, pass in references to external objects instead of temporaries:
 *
 * ```cpp
 * int i = 7;
 * auto e = foo | view::deep_take(i);
 * ```
 *
 * Wrapping an adaptor including its arguments:
 *
 * ```cpp
 *
 * std::vector<dna5_vector> foo{"AAATTT"_dna5, "CCCGGG"_dna5};
 *
 * auto t = foo | ranges::view::take(1);             // == [ [A,A,A,T,T,T] ]
 *
 * auto d = foo | view::deep{ranges::view::take(1)}; // == [ [A], [C] ]
 * // constructor arguments passed via {} and arguments to underlying view hardcoded inside
 *
 * // Or with an alias
 * namespace view
 * {
 * inline auto const deep_take1 = deep{ranges::view::take(1)};
 * }
 *
 * auto e = foo | view::deep_take1;                  // == [ [A], [C] ]
 * ```
 * In the above example the argument to the underlying adaptor is hardcoded and can't be changed at the call-site. It
 * is less flexible, but does not require workarounds for arguments that are expensive (or impossible) to copy.
 *
 */

template <typename underlying_adaptor_t>
class deep : public detail::pipable_adaptor_base<deep<underlying_adaptor_t>>
{
private:
    //!\brief Hold the inner functor by value or reference depending on deduced type.
    underlying_adaptor_t inner_functor;

public:

    /*!\name Constructors, destructor and assignment
     * \{
     */
    constexpr deep() = default;
    constexpr deep(deep const &) = default;
    constexpr deep(deep &&) = default;
    constexpr deep & operator=(deep const &) = default;
    constexpr deep & operator=(deep &&) = default;
    ~deep() = default;

    //!\brief Constructor that takes the underlying adaptor as parameter.
    constexpr deep(underlying_adaptor_t && inner)
        : inner_functor{std::forward<underlying_adaptor_t>(inner)}
    {}
    //!\}

    //!\privatesection

    /*!\brief       The implementation that the base classes operators all resolve to.
     * \tparam      urng_t Type of the underlying range.
     * \tparam      arg_types The argument types, not that **temporaries will be copied on recursion!**
     * \param[in]   urange The view's underlying range.
     * \param[in]   args Further arguments (optional).
     * \returns     A view with the inner adaptor applied on the innermost ranges.
     *
     * \details
     *
     * Recurses and calls ranges::view::transform if the underlying range is a range-of-ranges.
     */

    template <typename urng_t, typename ...arg_types>
    constexpr auto impl(urng_t && urange, arg_types && ...args) const
    {
        if constexpr (input_range_concept<urng_t> && input_range_concept<reference_t<urng_t>>)
        {
            return urange | ranges::view::transform(
                [argos = std::tuple<arg_types...>(std::forward<arg_types>(args)...), this] (auto && subrange)
                {
                    return impl_expand(std::forward<decltype(subrange)>(subrange),
                                    argos,
                                    std::make_index_sequence<std::tuple_size_v<decltype(argos)>>{});
                });
        }
        else
        {
            return inner_functor(std::forward<urng_t>(urange), std::forward<arg_types>(args)...);
        }
    }

    //!\brief Auxiliary function to expand formerly packed tuple.
    template <typename urng_t, typename tuple_t, size_t ... indexes>
    constexpr auto impl_expand(urng_t && urange, tuple_t && tuple, std::index_sequence<indexes...>) const
    {
        return impl(std::forward<urng_t>(urange), std::get<indexes>(std::forward<tuple_t>(tuple))...);
    }

};

/*!\name Template argument deduction guides.
 * \{
 */
//!\brief Template argument deduction helper that preserves lvalue references and turns rvalue references into values.
//!\relates deep
template <typename underlying_adaptor_t>
deep(underlying_adaptor_t && inner) -> deep<underlying_adaptor_t>;

//!\}

} // namespace seqan3::view
