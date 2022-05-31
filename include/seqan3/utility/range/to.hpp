// -----------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2022, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2022, Knut Reinert & MPI für molekulare Genetik
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at: https://github.com/seqan/seqan3/blob/master/LICENSE.md
// -----------------------------------------------------------------------------------------------------

/*!\file
 * \author Simon Gene Gottlieb <simon.gottlieb AT fu-berlin.de>
 * \brief Provides seqan3::ranges::to.
 */

#pragma once

#include <algorithm>
#include <ranges>

#include <seqan3/core/range/detail/adaptor_from_functor.hpp>
#include <seqan3/utility/concept.hpp>

namespace seqan3::detail
{

//!\brief Function object to convert a std::ranges::input_range to a fully defined container.
template <typename container_t>
    requires (!std::ranges::view<container_t>)
struct to_fn
{
private:
    /*!\brief Copies a range into a container.
     * \tparam rng_t       Type of the range.
     * \tparam container_t Type of the target container.
     */
    template <std::ranges::input_range rng_t>
        requires std::convertible_to<std::ranges::range_reference_t<rng_t>, std::ranges::range_value_t<container_t>>
    auto impl(rng_t && rng, container_t & container) const
    {
        std::ranges::copy(rng, std::back_inserter(container));
    }

    //!\overload
    template <std::ranges::input_range rng_t>
    auto impl(rng_t && rng, container_t & container) const
        requires std::ranges::input_range<std::ranges::range_value_t<rng_t>>
    {
        auto adaptor = to_fn<std::ranges::range_value_t<container_t>>{};
        auto inner_rng = rng | std::views::transform(adaptor);
        std::ranges::copy(inner_rng, std::back_inserter(container));
    }

public:
    /*!\brief Converts a template-template into a container.
     * \tparam rng_t  The type of the range being processed.
     * \tparam args_t The types of the arguments for the constructor.
     * \param  rng    The range being processed.
     * \param  args   Arguments to pass to the constructor of the container.
     */
    template <std::ranges::input_range rng_t, typename... args_t>
    constexpr auto operator()(rng_t && rng, args_t &&... args) const
    {
        auto new_container = container_t(std::forward<args_t>(args)...);

        // reserve memory if functionality is available
        if constexpr (std::ranges::sized_range<rng_t> && requires (container_t c) { c.reserve(size_t{}); })
            new_container.reserve(std::ranges::size(rng));

        impl(std::forward<rng_t>(rng), new_container);
        return new_container;
    }
};

/*!\brief Similar to to_fn, but accepts a template-template as argument,
 *        e.g.: to_fn<vector> instead of to_fn<vector<int>>.
 */
template <template <typename> typename container_t>
struct to_template_template_fn
{
    /*!\brief Converts a template-template into a container.
     * \tparam rng_t  The type of the range being processed.
     * \tparam args_t The types of the arguments for the constructor.
     * \param  rng    The range being processed.
     * \param  args   Arguments to pass to the constructor of the container.
     */
    template <std::ranges::input_range rng_t, typename... args_t>
    constexpr auto operator()(rng_t && rng, args_t &&... args) const
    {
        auto adaptor = to_fn<container_t<std::ranges::range_value_t<rng_t>>>{};
        return adaptor(std::forward<rng_t>(rng), std::forward<args_t>(args)...);
    }
};

} // namespace seqan3::detail

namespace seqan3::ranges
{

/*!\brief Converts a range to a container.
 * \ingroup utility_range
 * \details
 * To convert a range to a container, either the "pipe syntax" or the "function call" syntax can be used.
 * Both syntaxes support the explicit specification of the target container or
 * a specification with a deduced value type.
 *
 * Currently supported containers are:
 * - `std::vector`
 * - `std::list`
 * - `std::deque`
 *
 * \include doc/tutorial/05_ranges/to.cpp
 *
 * \noapi{This is a implementation of the C++23 ranges::to. It will be replaced with std::ranges::to.}
 */
template <typename container_t, typename... args_t>
constexpr auto to(args_t &&... args)
{
    return detail::adaptor_from_functor{detail::to_fn<container_t>{}, std::forward<args_t>(args)...};
}

/*!\brief Converts a range to a container.
 * \ingroup utility_range
 * \overload
 */
template <template <typename...> typename container_t, typename... args_t>
constexpr auto to(args_t &&... args)
{
    return detail::adaptor_from_functor{detail::to_template_template_fn<container_t>{}, std::forward<args_t>(args)...};
}

/*!\brief Converts a range to a container.
 * \ingroup utility_range
 * \overload
 */
template <typename container_t, std::ranges::input_range rng_t, typename... args_t>
constexpr auto to(rng_t && rng, args_t &&... args)
{
    return detail::adaptor_from_functor{detail::to_fn<container_t>{},
                                        std::forward<args_t>(args)...}(std::forward<rng_t>(rng));
}

/*!\brief Converts a range to a container.
 * \ingroup utility_range
 * \overload
 */
template <template <class...> typename container_t, std::ranges::input_range rng_t, typename... args_t>
constexpr auto to(rng_t && rng, args_t &&... args)
{
    return detail::adaptor_from_functor{detail::to_template_template_fn<container_t>{},
                                        std::forward<args_t>(args)...}(std::forward<rng_t>(rng));
}

} // namespace seqan3::ranges

namespace seqan3::views
{

/*!\brief Converts a range to a container.
 * \ingroup utility_views
 * \deprecated Use seqan3::ranges::to instead.
 */
template <typename container_t, typename... args_t>
SEQAN3_DEPRECATED_330 constexpr auto to(args_t &&... args)
{
    return ranges::to<container_t>(std::forward<args_t>(args)...);
}

/*!\brief Converts a range to a container.
 * \ingroup utility_views
 * \deprecated Use seqan3::ranges::to instead.
 */
template <template <typename...> typename container_t, typename... args_t>
SEQAN3_DEPRECATED_330 constexpr auto to(args_t &&... args)
{
    return ranges::to<container_t>(std::forward<args_t>(args)...);
}

/*!\brief Converts a range to a container.
 * \ingroup utility_views
 * \deprecated Use seqan3::ranges::to instead.
 */
template <typename container_t, std::ranges::input_range rng_t, typename... args_t>
SEQAN3_DEPRECATED_330 constexpr auto to(rng_t && rng, args_t &&... args)
{
    return ranges::to<container_t>(std::forward<rng_t>(rng), std::forward<args_t>(args)...);
}

/*!\brief Converts a range to a container.
 * \ingroup utility_views
 * \deprecated Use seqan3::ranges::to instead.
 */
template <template <typename...> typename container_t, std::ranges::input_range rng_t, typename... args_t>
SEQAN3_DEPRECATED_330 constexpr auto to(rng_t && rng, args_t &&... args)
{
    return ranges::to<container_t>(std::forward<rng_t>(rng), std::forward<args_t>(args)...);
}

} // namespace seqan3::views
