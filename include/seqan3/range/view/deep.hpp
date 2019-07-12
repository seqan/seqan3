// -----------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2019, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2019, Knut Reinert & MPI für molekulare Genetik
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at: https://github.com/seqan/seqan3/blob/master/LICENSE.md
// -----------------------------------------------------------------------------------------------------

/*!\file
 * \author Hannes Hauswedell <hannes.hauswedell AT fu-berlin.de>
 * \brief Provides seqan3::view::deep.
 */

#pragma once

#include <seqan3/range/view/detail.hpp>
#include <seqan3/std/ranges>

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
 * call std::view::reverse on a range-of-dna-ranges, it will revert *the order* of the dna-ranges, but leave
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
 * | range concepts and reference_t  | `urng_t` (underlying range type)      | `rrng_t` (returned range type)                     |
 * |---------------------------------|:-------------------------------------:|:--------------------------------------------------:|
 * | std::ranges::InputRange         | *required*                            | *preserved*                                        |
 * | std::ranges::ForwardRange       |                                       | *preserved*                                        |
 * | std::ranges::BidirectionalRange |                                       | *preserved*                                        |
 * | std::ranges::RandomAccessRange  |                                       | *preserved*                                        |
 * | std::ranges::ContiguousRange    |                                       | *lost*                                             |
 * |                                 |                                       |                                                    |
 * | std::ranges::ViewableRange      | *required*                            | *guaranteed*                                       |
 * | std::ranges::View               |                                       | *guaranteed*                                       |
 * | std::ranges::SizedRange         |                                       | *preserved*                                        |
 * | std::ranges::CommonRange        |                                       | *preserved*                                        |
 * | std::ranges::OutputRange        |                                       | *lost*                                             |
 * | seqan3::ConstIterableRange      |                                       | *preserved*                                        |
 * |                                 |                                       |                                                    |
 * | seqan3::reference_t             | std::ranges::InputRange           | std::ranges::InputRange + std::ranges::View |
 *
 * ### Examples
 *
 * Wrapping an adaptor that takes no parameters ("range adaptor <i>closure</i> object"):
 *
 * \snippet test/snippet/range/view/deep.cpp no_param
 *
 * Wrapping an adaptor that takes parameters:
 *
 * \snippet test/snippet/range/view/deep.cpp with_param
 *
 * The above example illustrates that view::deep has two sets of arguments, the **arguments to construct** this adaptor,
 * and the arguments passed to the underlying adaptor when **calling** this adaptor. You can use `()` for both, but
 * we highly recommend to use `{}` to not confuse these; or just use an alias.
 *
 * \attention Note that in the case of parameter handling the arguments to view::deep are **copied** to each invocation
 * of the underlying adaptor if they are temporaries. This is no problem for small objects like the integer above,
 * but might be expensive for larger ones. To avoid this, pass in references to external objects instead of temporaries:
 *
 * \snippet test/snippet/range/view/deep.cpp pass_ref
 *
 * Wrapping an adaptor including its arguments:
 *
 * \snippet test/snippet/range/view/deep.cpp wrap_args
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
    constexpr deep()                         noexcept = default; //!< Defaulted.
    constexpr deep(deep const &)             noexcept = default; //!< Defaulted.
    constexpr deep(deep &&)                  noexcept = default; //!< Defaulted.
    constexpr deep & operator=(deep const &) noexcept = default; //!< Defaulted.
    constexpr deep & operator=(deep &&)      noexcept = default; //!< Defaulted.
    ~deep()                                  noexcept = default; //!< Defaulted.

    using base_type::base_type;
    //!\}

    //!\privatesection

    //!\brief Explicitly inherit operator() because we are also defining our own.
    using base_type::operator();

    /*!\brief Unwrap the internal adaptor closure object and pipe the range into it.
     * \tparam    urng_t                Type of the underlying range.
     * \tparam    underlying_adaptor_t_ Same as underlying_adaptor_t with possibly different cref qualification.
     * \param[in] urange                The view's underlying range.
     * \param[in] adap                  The stored adaptor unwrapped by the base-classes explode()-member.
     * \returns A view with the inner adaptor applied on the innermost ranges.
     */
    template <std::ranges::InputRange urng_t, typename underlying_adaptor_t_>
    static constexpr auto impl(urng_t && urange, underlying_adaptor_t_ && adap)
    {
        return std::forward<urng_t>(urange) | std::forward<underlying_adaptor_t_>(adap);
    }

    /*!\brief Specialisation of the range-handling `operator()` for range-of-range (this is where deep
     * changes the behaviour for nested ranges).
     * \tparam    urng_t Type of the underlying range.
     * \param[in] urange The view's underlying range.
     * \returns A view with the inner adaptor applied on the innermost ranges.
     *
     * \details
     *
     * Recurses and calls std::view::transform if the underlying range is a range-of-ranges.
     */
    template <std::ranges::InputRange urng_t>
    //!\cond
        requires std::ranges::InputRange<reference_t<urng_t>>
    //!\endcond
    constexpr auto operator()(urng_t && urange) const &
    {
        return std::forward<urng_t>(urange) | std::view::transform([me = *this] (auto && e)
        {
            return std::forward<decltype(e)>(e) | me;
        });
    }

    //!\overload
    template <std::ranges::InputRange urng_t>
    //!\cond
        requires std::ranges::InputRange<reference_t<urng_t>>
    //!\endcond
    constexpr auto operator()(urng_t && urange) &&
    {
        return std::forward<urng_t>(urange) | std::view::transform([me = std::move(*this)] (auto && e)
        {
            return std::forward<decltype(e)>(e) | me;
        });
    }

    /*!\brief Called to produce a range adaptor closure object if the wrapped functor was **not** a range
     * adaptor *closure* object before, i.e. requires arguments.
     * \tparam    first_arg_t      The type of the first argument; must not model std::ranges::InputRange.
     * \tparam    stored_arg_types The argument types, note that **temporaries will be copied on recursion!**
     * \param[in] first            First argument.
     * \param[in] args             Further arguments (optional).
     * \returns A view::deep specialisation of a closure object, i.e. with stored arguments.
     */
    template <typename first_arg_t, typename ... stored_arg_types>
    //!\cond
        requires !std::ranges::InputRange<first_arg_t>
    //!\endcond
    constexpr auto operator()(first_arg_t && first, stored_arg_types && ...args) const
    {
        // The adaptor currently wrapped is a proto-adaptor and this function has the arguments to "complete" it.
        // We extract the adaptor that is stored and invoke it with the given arguments.
        // This returns an adaptor closure object.
        auto adaptor_closure = std::get<0>(this->arguments)(std::forward<first_arg_t>(first),
                                                            std::forward<stored_arg_types>(args)...);
        // Now we wrap this closure object back into a view::deep to get the deep behaviour.
        return deep<decltype(adaptor_closure)>{std::move(adaptor_closure)};
    }

    //!\overload
    constexpr auto operator()() const
    {
        // Proto-adaptors require arguments by definition, but some support defaulting those (e.g. view::translate).
        // This extracts the proto adaptor and invokes it without args which yields a different object, the closure
        // with defaulted arguments.
        auto adaptor_closure = std::get<0>(this->arguments)();
        // Now we wrap this closure object back into a view::deep to get the deep behaviour.
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
     * Recurses and calls std::view::transform if the underlying range is a range-of-ranges.
     */
    template <std::ranges::InputRange urng_t, typename ... stored_arg_types>
    //!\cond
        requires sizeof...(stored_arg_types) > 0
    //!\endcond
    constexpr auto operator()(urng_t && urange, stored_arg_types && ...args) const
    {
        auto adaptor_closure = std::get<0>(this->arguments)(std::forward<stored_arg_types>(args)...);
        return std::forward<urng_t>(urange) | std::move(adaptor_closure);
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
