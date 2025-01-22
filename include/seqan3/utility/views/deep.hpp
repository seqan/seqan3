// SPDX-FileCopyrightText: 2006-2025 Knut Reinert & Freie Universität Berlin
// SPDX-FileCopyrightText: 2016-2025 Knut Reinert & MPI für molekulare Genetik
// SPDX-License-Identifier: BSD-3-Clause

/*!\file
 * \brief Provides seqan3::views::deep.
 * \author Hannes Hauswedell <hannes.hauswedell AT fu-berlin.de>
 */

#pragma once

#include <ranges>

#include <seqan3/core/range/detail/adaptor_base.hpp>
#include <seqan3/core/range/type_traits.hpp>

namespace seqan3::views
{

/*!\brief   A wrapper type around an existing view adaptor that enables "deep view" behaviour for that view.
 * \ingroup utility_views
 * \tparam  underlying_adaptor_t The type of the adaptor being wrapped.
 *
 * \details
 *
 * ### Deep views
 *
 * If you pass a range to a view that view performs some transformation on that range. If the range passed is
 * multi-dimensional (i.e. a range-of-ranges) that transformation happens on the outermost range. So if you
 * call std::views::reverse on a range-of-dna-ranges, it will revert *the order* of the dna-ranges, but leave
 * the dna-ranges themselves unchanged.
 *
 * In some cases this is not desirable or even possible, i.e. seqan3::views::complement performs it's operation on
 * nucleotide-ranges and it would be logical to do so, even it is passed a range-of-nucleotide-ranges (it obviously
 * cannot transform the outer range). We call these views "deep views" as they always perform their operation on
 * the innermost ranges of a multi-dimensional range; in case the input is a one-dimensional range, deepness
 * does not modify the behaviour.
 *
 * ### Using views::deep
 *
 * Strictly speaking, seqan3::views::deep is a view adaptor adaptor, i.e. it gets passed **another adaptor when being
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
 * | Concepts and traits              | `urng_t` (underlying range type)      | `rrng_t` (returned range type)                     |
 * |----------------------------------|:-------------------------------------:|:--------------------------------------------------:|
 * | std::ranges::input_range         | *required*                            | *preserved*                                        |
 * | std::ranges::forward_range       |                                       | *preserved*                                        |
 * | std::ranges::bidirectional_range |                                       | *preserved*                                        |
 * | std::ranges::random_access_range |                                       | *preserved*                                        |
 * | std::ranges::contiguous_range    |                                       | *lost*                                             |
 * |                                  |                                       |                                                    |
 * | std::ranges::viewable_range      | *required*                            | *guaranteed*                                       |
 * | std::ranges::view                |                                       | *guaranteed*                                       |
 * | std::ranges::sized_range         |                                       | *preserved*                                        |
 * | std::ranges::common_range        |                                       | *preserved*                                        |
 * | std::ranges::output_range        |                                       | *lost*                                             |
 * | seqan3::const_iterable_range     |                                       | *preserved*                                        |
 * |                                  |                                       |                                                    |
 * | std::ranges::range_reference_t   | std::ranges::input_range              | std::ranges::input_range + std::ranges::view       |
 *
 * ### Examples
 *
 * Wrapping an adaptor that takes no parameters ("range adaptor <i>closure</i> object"):
 *
 * \include test/snippet/utility/views/deep_no_param.cpp
 *
 * Wrapping an adaptor that takes parameters:
 *
 * \include test/snippet/utility/views/deep_with_param.cpp
 *
 * The above example illustrates that views::deep has two sets of arguments, the **arguments to construct** this adaptor,
 * and the arguments passed to the underlying adaptor when **calling** this adaptor. You can use `()` for both, but
 * we highly recommend to use `{}` to not confuse these; or just use an alias.
 *
 * \attention Note that in the case of parameter handling the arguments to views::deep are **copied** to each invocation
 * of the underlying adaptor if they are temporaries. This is no problem for small objects like the integer above,
 * but might be expensive for larger ones. To avoid this, pass in references to external objects instead of temporaries:
 *
 * \include test/snippet/utility/views/deep_pass_ref.cpp
 *
 * Wrapping an adaptor including its arguments:
 *
 * \include test/snippet/utility/views/deep_wrap_args.cpp
 * In the above example the argument to the underlying adaptor is hardcoded and can't be changed at the call-site. It
 * is less flexible, but does not require workarounds for arguments that are expensive (or impossible) to copy.
 *
 */
template <typename underlying_adaptor_t>
class deep : public detail::adaptor_base<deep<underlying_adaptor_t>, underlying_adaptor_t>
{
private:
    //!\brief Type of the CRTP-base.
    using base_type = detail::adaptor_base<deep<underlying_adaptor_t>, underlying_adaptor_t>;

    //!\brief Befriend the base class so it can call impl().
    friend base_type;

public:
    /*!\name Constructors, destructor and assignment
     * \{
     */
    constexpr deep() noexcept = default;                         //!< Defaulted.
    constexpr deep(deep const &) noexcept = default;             //!< Defaulted.
    constexpr deep(deep &&) noexcept = default;                  //!< Defaulted.
    constexpr deep & operator=(deep const &) noexcept = default; //!< Defaulted.
    constexpr deep & operator=(deep &&) noexcept = default;      //!< Defaulted.
    ~deep() noexcept = default;                                  //!< Defaulted.

    using base_type::base_type;
    //!\}

    //!\privatesection

    //!\brief Explicitly inherit operator() because we are also defining our own.
    using base_type::operator();

    /*!\brief Unwrap the internal adaptor closure object and pipe the range into it.
     * \tparam    urng_t                Type of the underlying range.
     * \tparam    underlying_adaptor_t_ Same as underlying_adaptor_t with possibly different cref qualification.
     * \param[in] urange                The view's underlying range.
     * \param[in] deep_adaptor          The stored adaptor unwrapped by the base-classes explode()-member.
     * \returns A view with the inner adaptor applied on the innermost ranges.
     */
    template <std::ranges::input_range urng_t, typename underlying_adaptor_t_>
    static constexpr auto impl(urng_t && urange, underlying_adaptor_t_ && deep_adaptor)
    {
        static_assert(std::same_as<underlying_adaptor_t_, underlying_adaptor_t>,
                      "Internally stored deep-adaptor does not match!");

        constexpr size_t range_dimension = range_dimension_v<urng_t>;

        // note: if range_dimension == 1, the expression will actually be
        // `std::forward<underlying_adaptor_t_>(adap)(urange)`, thus allowing the stateful adaptor to move its arguments
        // into the view that will be constructed.
        return recursive_adaptor<range_dimension>(std::forward<underlying_adaptor_t_>(deep_adaptor))(urange);
    }

    /*!\brief recursively construct the deep adaptor
     * \details
     *
     * this function recursively constructs a range adaptor:
     *
     * ```cpp
     * return std::views::transform([](auto && inner_range_depth1)
     * {
     *     return inner_range_depth1 | std::views::transform([](auto && inner_range_depth2)
     *     {
     *         // ...
     *         return inner_range_depth_n | std::views::transform([adap](auto && innermost_range)
     *         {
     *             // only innermost std::views::transform stores the actual deep-adaptor (stored only once!)
     *             return innermost_range | std::forward<underlying_adaptor_t>(deep_adaptor);
     *         });
     *         // ...
     *     });
     * });
     * ```
     */
    template <std::size_t range_dimension>
    static constexpr decltype(auto) recursive_adaptor(underlying_adaptor_t deep_adaptor)
    {
        if constexpr (range_dimension > 1u)
        {
            auto transform = [adaptor = recursive_adaptor<range_dimension - 1u>(
                                  std::forward<underlying_adaptor_t>(deep_adaptor))](auto && inner_range)
            {
                // We don't want to move a stateful adaptor here, as this adaptor will be called on any element of this
                // std::views::transform range.
                return adaptor(std::forward<decltype(inner_range)>(inner_range));
            };
            return std::views::transform(std::move(transform));
        }
        else
        {
            // recursion anchor: only the innermost std::views::transform will store the deep_adaptor.
            // In an earlier version of seqan3 we recursively passed the deep_adaptor through all std::view::transform
            // instances and each depth stored the same deep_adaptor. The current approach can save memory as it only
            // stores the deep_adaptor once in the innermost std::views::transform. It will also produce slightly better
            // stack-traces by defining the recursion over the range_dimension.

            // NOTE: to allow an lvalue to be returned we need to have the return type decltype(auto).
            return deep_adaptor;
        }
    }

    /*!\brief Called to produce a range adaptor closure object if the wrapped functor was **not** a range
     * adaptor *closure* object before, i.e. requires arguments.
     * \tparam    first_arg_t      The type of the first argument; must not model std::ranges::input_range.
     * \tparam    stored_arg_types The argument types, note that **temporaries will be copied on recursion!**
     * \param[in] first            First argument.
     * \param[in] args             Further arguments (optional).
     * \returns A views::deep specialisation of a closure object, i.e. with stored arguments.
     */
    template <typename first_arg_t, typename... stored_arg_types>
        requires (!std::ranges::input_range<first_arg_t>)
    constexpr auto operator()(first_arg_t && first, stored_arg_types &&... args) const
    {
        // The adaptor currently wrapped is a proto-adaptor and this function has the arguments to "complete" it.
        // We extract the adaptor that is stored and invoke it with the given arguments.
        // This returns an adaptor closure object.
        auto adaptor_closure =
            std::get<0>(this->arguments)(std::forward<first_arg_t>(first), std::forward<stored_arg_types>(args)...);
        // Now we wrap this closure object back into a views::deep to get the deep behaviour.
        return deep<decltype(adaptor_closure)>{std::move(adaptor_closure)};
    }

    //!\overload
    constexpr auto operator()() const
    {
        // Proto-adaptors require arguments by definition, but some support defaulting those (e.g. views::translate).
        // This extracts the proto adaptor and invokes it without args which yields a different object, the closure
        // with default arguments.
        auto adaptor_closure = std::get<0>(this->arguments)();
        // Now we wrap this closure object back into a views::deep to get the deep behaviour.
        return deep<decltype(adaptor_closure)>{std::move(adaptor_closure)};
    }

    /*!\brief Called to produce a range if the wrapped functor was **not** a range adaptor *closure* object
     * before but necessary arguments are also provided.
     * \tparam    urng_t    Type of the underlying range.
     * \tparam    arg_types The argument types, note that **temporaries will be copied on recursion!**
     * \param[in] urange    The view's underlying range.
     * \param[in] args      Further arguments.
     * \returns A view with the inner adaptor applied on the innermost ranges.
     *
     * \details
     *
     * Recurses and calls std::views::transform if the underlying range is a range-of-ranges.
     */
    template <std::ranges::input_range urng_t, typename... stored_arg_types>
        requires (sizeof...(stored_arg_types) > 0)
    constexpr auto operator()(urng_t && urange, stored_arg_types &&... args) const
    {
        auto adaptor_closure = std::get<0>(this->arguments)(std::forward<stored_arg_types>(args)...);
        deep<decltype(adaptor_closure)> deep_adaptor{std::move(adaptor_closure)};
        return deep_adaptor(std::forward<urng_t>(urange));
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

} // namespace seqan3::views
