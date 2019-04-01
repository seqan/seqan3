// -----------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2019, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2019, Knut Reinert & MPI für molekulare Genetik
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at: https://github.com/seqan/seqan3/blob/master/LICENSE
// -----------------------------------------------------------------------------------------------------

/*!\file
 * \author Hannes Hauswedell <hannes.hauswedell AT fu-berlin.de>
 * \brief Provides seqan3::view::drop.
 */

#pragma once

#include <seqan3/io/exception.hpp>
#include <seqan3/range/concept.hpp>
#include <seqan3/range/view/detail.hpp>
#include <seqan3/std/concepts>
#include <seqan3/std/iterator>
#include <seqan3/std/ranges>
#include <seqan3/std/span>
#include <seqan3/std/type_traits>

namespace seqan3::detail
{

// ============================================================================
//  drop_fn (adaptor definition)
// ============================================================================

/*!\brief View adaptor definition for view::drop and view::drop_or_throw.
 */
class drop_fn : public pipable_adaptor_base<drop_fn>
{
private:
    //!\brief Type of the CRTP-base.
    using base_t = pipable_adaptor_base<drop_fn>;

public:
    //!\brief Inherit the base class's Constructors.
    using base_t::base_t;

private:
    //!\brief Befriend the base class so it can call impl().
    friend base_t;

    /*!\brief       Call the view's constructor with the underlying view as argument.
     * \returns     An instance of std::ranges::subrange.
     */
    template <std::ranges::ViewableRange urng_t>
    static auto impl(urng_t && urange, size_t drop_size)
    {
        auto b = std::ranges::begin(urange);

        if constexpr (std::ranges::SizedRange<urng_t>)
        {
            drop_size = std::min(drop_size, std::ranges::size(urange));
            std::ranges::advance(b, drop_size);

            return std::ranges::subrange<std::ranges::iterator_t<urng_t>,
                                         std::ranges::sentinel_t<urng_t>,
                                         std::ranges::subrange_kind::sized>
            {
                std::move(b),
                std::ranges::end(urange),
                std::ranges::size(urange) - drop_size
            };
        }
        else
        {
            for (size_t i = 0; (i < drop_size) && (b != std::ranges::end(urange)); ++i, ++b);

            return std::ranges::subrange<std::ranges::iterator_t<urng_t>,
                                         std::ranges::sentinel_t<urng_t>>
            {
                std::move(b),
                std::ranges::end(urange)
            };
        }
    }

    /*!\brief       Overload for contiguous, sized ranges.
     * \returns     A std::span over the input.
     */
    template <std::ranges::ViewableRange urng_t>
    //!\cond
        requires std::ranges::ContiguousRange<urng_t> && std::ranges::SizedRange<urng_t>
    //!\endcond
    static auto impl(urng_t && urange, size_t drop_size)
    {
        drop_size = std::min(drop_size, std::ranges::size(urange));
        return std::span{std::ranges::data(urange) + drop_size, std::ranges::size(urange) - drop_size};
    }

    /*!\brief       Overload for std::basic_string_view.
     * \returns     A std::basic_string_view over the input.
     */
    template <std::ranges::ViewableRange urng_t>
    //!\cond
        requires std::ranges::ContiguousRange<urng_t> && std::ranges::SizedRange<urng_t> &&
                 is_type_specialisation_of_v<std::remove_reference_t<urng_t>, std::basic_string_view>
    //!\endcond
    static auto impl(urng_t && urange, size_t drop_size)
    {
        drop_size = std::min(drop_size, std::ranges::size(urange));
        return urange.substr(drop_size);
    }

    /*!\brief       Overload for std::basic_string.
     * \returns     A std::basic_string_view over the input.
     */
    template <typename char_t, typename traits_t, typename alloc_t>
    static std::basic_string_view<char_t, traits_t> impl(std::basic_string<char_t, traits_t, alloc_t> & urange,
                                                         size_t drop_size)
    {
        drop_size = std::min(drop_size, std::ranges::size(urange));
        return {std::ranges::data(urange) + drop_size, std::ranges::size(urange) - drop_size};
    }

    //!\overload
    template <typename char_t, typename traits_t, typename alloc_t>
    static std::basic_string_view<char_t, traits_t> impl(std::basic_string<char_t, traits_t, alloc_t> const & urange,
                                                         size_t drop_size)
    {
        drop_size = std::min(drop_size, std::ranges::size(urange));
        return {std::ranges::data(urange) + drop_size, std::ranges::size(urange) - drop_size};
    }

    //!\overload
    template <typename char_t, typename traits_t, typename alloc_t>
    static std::basic_string_view<char_t, traits_t> impl(std::basic_string<char_t, traits_t, alloc_t> const && urange,
                                                         size_t drop_size) = delete;
};

} // namespace seqan3::detail

// ============================================================================
//  view::drop (adaptor instance definition)
// ============================================================================

namespace seqan3::view
{

/*!\name General purpose views
 * \{
 */

/*!\brief               A view adaptor that returns all elements after n from the underlying range (or an empty range
 *                      if the underlying range is shorter).
 * \tparam urng_t       The type of the range being processed. See below for requirements. [template parameter is
 *                      omitted in pipe notation]
 * \param[in] urange    The range being processed. [parameter is omitted in pipe notation]
 * \param[in] drop_size The number of elements to drop from the beginning.
 * \returns             All elements of the underlying range after the first `drop_size`.
 * \ingroup view
 *
 * \details
 *
 * ### View properties
 *
 * | range concepts and reference_t  | `urng_t` (underlying range type)      | `rrng_t` (returned range type)                     |
 * |---------------------------------|:-------------------------------------:|:--------------------------------------------------:|
 * | std::ranges::InputRange         | *required*                            | *preserved*                                        |
 * | std::ranges::ForwardRange       |                                       | *preserved*                                        |
 * | std::ranges::BidirectionalRange |                                       | *preserved*                                        |
 * | std::ranges::RandomAccessRange  |                                       | *preserved*                                        |
 * | std::ranges::ContiguousRange    |                                       | *preserved*                                        |
 * |                                 |                                       |                                                    |
 * | std::ranges::ViewableRange      | *required*                            | *guaranteed*                                       |
 * | std::ranges::View               |                                       | *guaranteed*                                       |
 * | std::ranges::SizedRange         |                                       | *preserved*                                        |
 * | std::ranges::CommonRange        |                                       | *preserved*                                        |
 * | std::ranges::OutputRange        |                                       | *preserved* unless `urng_t` is std::basic_string   |
 * | seqan3::const_iterable_concept  |                                       | *preserved*                                        |
 * |                                 |                                       |                                                    |
 * | seqan3::reference_t             |                                       | seqan3::reference_t<urng_t>                        |
 *
 * See the \link view view submodule documentation \endlink for detailed descriptions of the view properties.
 *
 * ### Return type
 *
 * | `urng_t` (underlying range type)                          | `rrng_t` (returned range type)                     |
 * |:---------------------------------------------------------:|:--------------------------------------------------:|
 * | std::basic_string *or* std::basic_string_view             | std::basic_string_view                             |
 * | std::ranges::SizedRange && std::ranges::ContiguousRange   | std::span                                          |
 * | *else*                                                    | std::ranges::subrange                              |
 *
 * The adaptor is different from std::view::drop in that it performs type erasure for some underlying ranges.
 * It returns exactly the type specified above.
 *
 * ### Complexity
 *
 * Construction time of the returned view is in \f$ O(1) \f$ if the underlying range models at least
 * std::ranges::RandomAccessRange and std::ranges::SizedRange; otherwise in \f$ O(drop\_size) \f$.
 *
 * ### Example
 *
 * \include test/snippet/range/view/drop.cpp
 *
 * \hideinitializer
 */
inline constexpr auto drop = detail::drop_fn{};

//!\}

} // namespace seqan3::view
