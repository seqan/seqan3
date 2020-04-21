// -----------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2020, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2020, Knut Reinert & MPI für molekulare Genetik
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at: https://github.com/seqan/seqan3/blob/master/LICENSE.md
// -----------------------------------------------------------------------------------------------------

/*!\file
 * \author Enrico Seiler <enrico.seiler AT fu-berlin.de>
 * \brief Provides seqan3::views::kmer_hash.
 */

#pragma once

#include <range/v3/view/sliding.hpp>

#include <seqan3/alphabet/concept.hpp>
#include <seqan3/core/type_traits/range.hpp>
#include <seqan3/range/views/detail.hpp>
#include <seqan3/std/ranges>

namespace seqan3::detail
{
// ============================================================================
//  naive_kmer_hash_fn (adaptor definition)
// ============================================================================

//!\brief views::kmer_hash's range adaptor object type (non-closure).
struct naive_kmer_hash_fn
{
    //!\brief Store the argument and return a range adaptor closure object.
    constexpr auto operator()(size_t const k) const noexcept
    {
        return detail::adaptor_from_functor{*this, k};
    }

    /*!\brief            Call the view's constructor with the underlying view as argument.
     * \param[in] urange The input range to process. Must model std::ranges::viewable_range and the reference type of the
     *                   range of the range must model seqan3::semialphabet.
     * \param[in] k      The k-mer size to construct hashes for.
     * \returns          A range of converted elements.
     */
    template <std::ranges::viewable_range urng_t>
    //!\cond
        requires semialphabet<std::ranges::range_reference_t<urng_t>>
    //!\endcond
    constexpr auto operator()(urng_t && urange, size_t const k) const noexcept
    {
        return std::forward<urng_t>(urange) | ranges::view::sliding(k) | std::views::transform(
        [] (auto const in)
        {
            std::hash<decltype(in)> h{};
            return h(in);
        });
    }
};

} // namespace seqan3::detail

namespace seqan3::views
{

/*!\name Alphabet related views
 * \{
 */

/*!\brief               A view that calls std::hash on each substring of length k in the input range.
 * \tparam urng_t       The type of the range being processed. See below for requirements. [template parameter is
 *                      omitted in pipe notation]
 * \param[in] urange    The range being processed. [parameter is omitted in pipe notation]
 * \returns             A range of unsigned integral values where each value is the hash of the resp. k-mer.
 *                      See below for the properties of the returned range.
 * \ingroup views
 *
 * ### View properties
 *
 * | Concepts and traits              | `urng_t` (underlying range type)      | `rrng_t` (returned range type)                     |
 * |----------------------------------|:-------------------------------------:|:--------------------------------------------------:|
 * | std::ranges::input_range         | *required*                            | *preserved*                                        |
 * | std::ranges::forward_range       | *required*                            | *preserved*                                        |
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
 * | std::ranges::range_reference_t   | seqan3::semialphabet                  | std::size_t                                        |
 *
 * See the \link views views submodule documentation \endlink for detailed descriptions of the view properties.
 *
 * ### Example
 * \snippet test/snippet/range/views/kmer_hash.cpp usage
 * \hideinitializer
 */
inline auto constexpr naive_kmer_hash = detail::naive_kmer_hash_fn{};

} // namespace seqan3::views
