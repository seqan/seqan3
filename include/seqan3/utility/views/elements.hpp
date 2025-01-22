// SPDX-FileCopyrightText: 2006-2025 Knut Reinert & Freie Universität Berlin
// SPDX-FileCopyrightText: 2016-2025 Knut Reinert & MPI für molekulare Genetik
// SPDX-License-Identifier: BSD-3-Clause

/*!\file
 * \brief Provides seqan3::views::elements.
 * \author Joshua Kim <joshua.kim AT fu-berlin.de>
 */

#pragma once

#include <ranges>

#include <seqan3/utility/tuple/concept.hpp>
#include <seqan3/utility/type_traits/basic.hpp>

namespace seqan3::views
{
/*!\brief               A view calling `get` on each element in a range.
 * \tparam size_t       The type of the range being processed. See below for requirements. [template parameter is
 *                      omitted in pipe notation]
 * \tparam index        The index to get.
 * \param[in] urange    The range being processed. [parameter is omitted in pipe notation]
 * \returns             A range of elements where every element is the result of calling `get<index>` on the underlying
 *                      element.
                        See below for the properties of the returned range.
 * \ingroup utility_views
 * \sa https://en.cppreference.com/w/cpp/ranges/elements_view
 *
 * \details
 *
 * \header_file{seqan3/utility/views/elements.hpp}
 *
 * This view may be used instead of `std::view::elements` when the underlying range contains seqan3 types, e.g.
 * seqan3::qualified.
 *
 * \if DEV
 * `std::views::elements` uses, per standard, a qualified call to `get`, i.e. `std::get` (see
 * [GCC ticket](https://gcc.gnu.org/bugzilla/show_bug.cgi?id=100233)).
 *
 * Combined with the way ADL does (not) work in this case, it is not possible to overload `std::get` for our types
 * ([Godbolt](https://godbolt.org/z/dGeqzxarv)).
 *
 * This is why we need a `seqan3::views::elements`; bear in mind that our requirement of `tuple_like` on the incoming
 * elements might not be completely equivalent to the STL's requirements
 * ([has-tuple-element](https://eel.is/c++draft/range.elements.view#concept:has-tuple-element)).
 * \endif
 * ### View properties
 *
 * | Concepts and traits              | `urng_t` (underlying range type)      | `rrng_t` (returned range type)                                     |
 * |----------------------------------|:-------------------------------------:|:------------------------------------------------------------------:|
 * | std::ranges::input_range         | *required*                            | *preserved*                                                        |
 * | std::ranges::forward_range       |                                       | *preserved*                                                        |
 * | std::ranges::bidirectional_range |                                       | *preserved*                                                        |
 * | std::ranges::random_access_range |                                       | *preserved*                                                        |
 * | std::ranges::contiguous_range    |                                       | *lost*                                                             |
 * |                                  |                                       |                                                                    |
 * | std::ranges::viewable_range      | *required*                            | *preserved*                                                        |
 * | std::ranges::view                |                                       | *preserved*                                                        |
 * | std::ranges::sized_range         |                                       | *preserved*                                                        |
 * | std::ranges::common_range        |                                       | *preserved*                                                        |
 * | std::ranges::output_range        |                                       | *preserved*                                                        |
 * | seqan3::const_iterable_range     |                                       | *preserved*                                                        |
 * |                                  |                                       |                                                                    |
 * | std::ranges::range_reference_t   | seqan3::tuple_like                    | std::tuple_element_t<index, std::ranges::range_reference_t<urng_t>>|
 *
 * See the \link views views submodule documentation \endlink for detailed descriptions of the view properties.
 *
 * ### Example
 *
 * \include test/snippet/utility/views/elements.cpp
 * \hideinitializer
 *
 * \experimentalapi{Experimental since version 3.1.}
 */
template <auto index>
inline constexpr auto elements = std::views::transform(
    [](auto && in) -> decltype(auto)
    {
        using std::get;
        using seqan3::get;

        using element_t = decltype(in);

        static_assert(tuple_like<element_t>,
                      "You may only pass ranges to views::element_t whose reference_t models tuple_like.");

        // we need to explicitly remove && around temporaries to return values as values (and not as rvalue references)
        // we cannot simply cast to std::tuple_element_t (or set that as return value), because some tuples, like
        // our alphabet_tuple_base alphabets do not return that type when get is called on them (they return a proxy)
        using ret_type = remove_rvalue_reference_t<decltype(get<index>(std::forward<element_t>(in)))>;
        return static_cast<ret_type>(get<index>(std::forward<element_t>(in)));
    });

} // namespace seqan3::views
