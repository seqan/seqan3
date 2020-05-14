// -----------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2020, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2020, Knut Reinert & MPI für molekulare Genetik
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at: https://github.com/seqan/seqan3/blob/master/LICENSE.md
// -----------------------------------------------------------------------------------------------------

/*!\file
 * \author Mitra Darvish <mitra.darvish AT fu-berlin.de>
 * \brief Provides seqan3::views::naive_minimiser_hash.
 */

#pragma once

#include <seqan3/std/ranges>

#include <range/v3/view/sliding.hpp>

#include <seqan3/alphabet/concept.hpp>
#include <seqan3/core/type_traits/range.hpp>
#include <seqan3/range/views/detail.hpp>
#include <seqan3/range/views/kmer_hash.hpp>
#include <seqan3/range/views/minimiser.hpp>

namespace seqan3::detail
{
// ============================================================================
//  naive_minimiser_hash_fn (adaptor definition)
// ============================================================================

//!\brief views::naive_minimiser_hash's range adaptor object type (non-closure).
struct naive_minimiser_hash_fn
{
    /*!\brief                Store the shape and the window size and return a range adaptor closure object.
    * \param[in] shape       The seqan3::shape to use for hashing.
    * \param[in] window_size The windows size to use.
    * \throws std::invalid_argument if the size of the shape is greater than the `window_size`.
    * \returns               A range of converted elements.
    */
    constexpr auto operator()(shape const & shape, uint32_t const window_size) const
    {
        return seqan3::detail::adaptor_from_functor{*this, shape, window_size};
    }

    /*!\brief                 Call the view's constructor with the underlying view, a seqan3::shape and a window size as
     *                        argument.
     * \param[in] urange      The input range to process. Must model std::ranges::viewable_range and the reference type
     *                        of the range must model seqan3::semialphabet.
     * \param[in] shape       The seqan3::shape to use for hashing.
     * \param[in] window_size The size of the window.
     * \param[in] seed        The seed to use.
     * \throws std::invalid_argument if the size of the shape is greater than the `window_size`.
     * \returns               A range of converted elements.
     */
   template <std::ranges::range urng_t>
   constexpr auto operator()(urng_t && urange, shape const & shape, uint32_t const window_size,
                             uint64_t const seed = 0x8F3F73B5CF1C9ADE) const
   {
       static_assert(std::ranges::viewable_range<urng_t>,
           "The range parameter to views::minimiser_hash cannot be a temporary of a non-view range.");
       static_assert(std::ranges::forward_range<urng_t>,
           "The range parameter to views::minimiser_hash must model std::ranges::forward_range.");
       static_assert(semialphabet<reference_t<urng_t>>,
           "The range parameter to views::minimiser_hash must be over elements of seqan3::semialphabet.");
       if (shape.size() > window_size)
           throw std::invalid_argument{"The size of the shape cannot be greater than the window size."};

       return std::forward<urng_t>(urange) | seqan3::views::kmer_hash(shape)
                                           | std::views::transform([seed] (uint64_t i) { return i ^ seed; })
                                           | ranges::view::sliding(window_size - shape.size() + 1)
                                           | std::views::transform([] (auto const in)
                                                                   {
                                                                       return *std::min_element(in.begin(), in.end());
                                                                   });
   }
};

} // namespace seqan3::detail

namespace seqan3::views
{

/*!\name Alphabet related views
 * \{
 */

/*!\brief               A view that calls std::min_element on each substring of length window_size - shape.size() + 1 in
 *                      the input range.
 * \tparam urng_t       The type of the range being processed. See below for requirements. [template parameter is
 *                      omitted in pipe notation]
 * \param[in] urange    The range being processed. [parameter is omitted in pipe notation]
 * \returns             A range of unsigned integral values where each value is the minimiser of the resp. window.
 *                      See below for the properties of the returned range.
 * \ingroup views
 *
 * ### View properties
 *
 * | Concepts and traits              | `urng_t` (underlying range type)      | `rrng_t` (returned range type)                     |
 * |----------------------------------|:-------------------------------------:|:--------------------------------------------------:|
 * | std::ranges::input_range         | *required*                            | *preserved*                                        |
 * | std::ranges::forward_range       | *required*                            | *preserved*                                        |
 * | std::ranges::bidirectional_range |                                       | *lost*                                             |
 * | std::ranges::random_access_range |                                       | *lost*                                             |
 * | std::ranges::contiguous_range    |                                       | *lost*                                             |
 * |                                  |                                       |                                                    |
 * | std::ranges::viewable_range      | *required*                            | *guaranteed*                                       |
 * | std::ranges::view                |                                       | *guaranteed*                                       |
 * | std::ranges::sized_range         |                                       | *lost*                                             |
 * | std::ranges::common_range        |                                       | *lost*                                             |
 * | std::ranges::output_range        |                                       | *lost*                                             |
 * | seqan3::const_iterable_range     |                                       | *preserved*                                        |
 * |                                  |                                       |                                                    |
 * | std::ranges::range_reference_t   | seqan3::semialphabet                  | std::size_t                                        |
 *
 * See the \link views submodule documentation \endlink for detailed descriptions of the view properties.
 *
 */
inline auto constexpr naive_minimiser_hash = detail::naive_minimiser_hash_fn{};

//!\}

} // namespace seqan3::views
