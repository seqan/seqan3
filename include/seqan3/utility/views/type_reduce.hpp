// SPDX-FileCopyrightText: 2006-2025 Knut Reinert & Freie Universität Berlin
// SPDX-FileCopyrightText: 2016-2025 Knut Reinert & MPI für molekulare Genetik
// SPDX-License-Identifier: BSD-3-Clause

/*!\file
 * \author Hannes Hauswedell <hannes.hauswedell AT fu-berlin.de>
 * \brief Provides seqan3::views::type_reduce.
 */

#pragma once

#include <concepts>
#include <ranges>
#include <span>
#include <string_view>

#include <seqan3/core/detail/template_inspection.hpp>
#include <seqan3/core/range/detail/adaptor_base.hpp>

namespace seqan3::detail
{

// ============================================================================
//  type_reduce_fn (adaptor definition)
// ============================================================================

/*!\brief View adaptor definition for views::type_reduce.
 */
class type_reduce_fn : public adaptor_base<type_reduce_fn>
{
private:
    //!\brief Type of the CRTP-base.
    using base_t = adaptor_base<type_reduce_fn>;

public:
    //!\brief Inherit the base class's Constructors.
    using base_t::base_t;

private:
    //!\brief Befriend the base class so it can call impl().
    friend base_t;

    /*!\brief Type erase if possible and delegate to std::views::all otherwise.
     * \returns An instance of std::span, std::basic_string_view, std::ranges::subrange or std::views::all's return type.
     */
    template <std::ranges::range urng_t>
    static constexpr auto impl(urng_t && urange)
    {
        static_assert(
            std::ranges::viewable_range<urng_t>,
            "The views::type_reduce adaptor can only be passed viewable_ranges, i.e. Views or &-to-non-View.");

        // string_view
        if constexpr (is_type_specialisation_of_v<std::remove_cvref_t<urng_t>, std::basic_string_view>)
        {
            return std::views::all(std::forward<urng_t>(urange));
        }
        // string const &
        else if constexpr (is_type_specialisation_of_v<std::remove_cvref_t<urng_t>, std::basic_string>
                           && std::is_const_v<std::remove_reference_t<urng_t>>)
        {
            return std::basic_string_view{std::ranges::data(urange), std::ranges::size(urange)};
        }
        // contiguous
        else if constexpr (std::ranges::borrowed_range<urng_t> && std::ranges::contiguous_range<urng_t>
                           && std::ranges::sized_range<urng_t>)
        {
            return std::span{std::ranges::data(urange), std::ranges::size(urange)};
        }
        // random_access
        else if constexpr (std::ranges::borrowed_range<urng_t> && std::ranges::random_access_range<urng_t>
                           && std::ranges::sized_range<urng_t>)
        {
            return std::ranges::subrange<std::ranges::iterator_t<urng_t>, std::ranges::iterator_t<urng_t>>{
                std::ranges::begin(urange),
                std::ranges::begin(urange) + std::ranges::size(urange),
                std::ranges::size(urange)};
        }
        // pass to std::views::all (will return a view or ref-view)
        else
        {
            return std::views::all(std::forward<urng_t>(urange));
        }
    }
};

} // namespace seqan3::detail

// ============================================================================
//  views::type_reduce (adaptor instance definition)
// ============================================================================

namespace seqan3::views
{
/*!\brief               A view adaptor that behaves like std::views::all, but type erases certain ranges.
 * \tparam urng_t       The type of the range being processed. See below for requirements. [template parameter is
 *                      omitted in pipe notation]
 * \param[in] urange    The range being processed. [parameter is omitted in pipe notation]
 * \returns             The range turned into a view.
 * \ingroup utility_views
 *
 * \details
 *
 * \header_file{seqan3/utility/views/view_type_reduce.hpp}
 *
 * ### View properties
 *
 * | Concepts and traits              | `urng_t` (underlying range type)  | `rrng_t` (returned range type)         |
 * |----------------------------------|:---------------------------------:|:--------------------------------------:|
 * | std::ranges::input_range         | *required*                        | *preserved*                            |
 * | std::ranges::forward_range       |                                   | *preserved*                            |
 * | std::ranges::bidirectional_range |                                   | *preserved*                            |
 * | std::ranges::random_access_range |                                   | *preserved*                            |
 * | std::ranges::contiguous_range    |                                   | *preserved*                            |
 * |                                  |                                   |                                        |
 * | std::ranges::viewable_range      | *required*                        | *guaranteed*                           |
 * | std::ranges::view                |                                   | *guaranteed*                           |
 * | std::ranges::sized_range         |                                   | *preserved*                            |
 * | std::ranges::common_range        |                                   | *preserved*                            |
 * | std::ranges::output_range        |                                   | *preserved*                            |
 * | seqan3::const_iterable_range     |                                   | *preserved*                            |
 * |                                  |                                   |                                        |
 * | std::ranges::range_reference_t   |                                   | std::ranges::range_reference_t<urng_t> |
 *
 * See the \link views views submodule documentation \endlink for detailed descriptions of the view properties.
 *
 * ### Return type
 *
 * | `urng_t` (underlying range type)                                                                  | `rrng_t` (returned range type)  |
 * |:-------------------------------------------------------------------------------------------------:|:-------------------------------:|
 * | `std::ranges::view`                                                                               | `urng_t`                        |
 * | `std::basic_string const &` *or* `std::basic_string_view`                                         | `std::basic_string_view`        |
 * | `std::ranges::borrowed_range && std::ranges::sized_range && std::ranges::contiguous_range`        | `std::span`                     |
 * | `std::ranges::borrowed_range && std::ranges::sized_range && std::ranges::random_access_range`     | `std::ranges::subrange`         |
 * | *else*                                                                                            | `std::ranges::ref_view<urng_t>` |
 *
 * This adaptor is different from std::views::all in that it performs type erasure for some underlying ranges;
 * std::views::all always returns `std::ranges::ref_view<urng_t>` or the type itself if it already is a view.
 *
 * ### Example
 *
* \include test/snippet/utility/views/type_reduce.cpp
 *
 * \hideinitializer
 * \experimentalapi{Experimental since version 3.1.}
 */
inline constexpr auto type_reduce = detail::type_reduce_fn{};

} // namespace seqan3::views

namespace seqan3
{
/*!\brief Deduces the return value of seqan3::views::type_reduce.
 * \details
 * \noapi
 * \sa https://en.cppreference.com/w/cpp/ranges/all_view
 */
template <typename t>
using type_reduce_t = decltype(views::type_reduce(std::declval<t>()));
} // namespace seqan3
