// -----------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2020, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2020, Knut Reinert & MPI für molekulare Genetik
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at: https://github.com/seqan/seqan3/blob/master/LICENSE.md
// -----------------------------------------------------------------------------------------------------

/*!\file
 * \author Mitra Darvish <mitra.darvish AT fu-berlin.de>
 * \brief Provides seqan3::views::minimizer_hash.
 */

#pragma once

#include <seqan3/range/views/kmer_hash.hpp>
#include <seqan3/range/views/minimizer.hpp>
#include <seqan3/search/kmer_index/shape.hpp>

namespace seqan3::detail
{
// ============================================================================
//  views::minimizer_hash (adaptor instance definition)
// ============================================================================

//![adaptor_def]
//!\brief views::minimizer's range adaptor object type (non-closure).
struct minimizer_hash_fn
{
    //!\brief Store the shape and the window_size and return a range adaptor closure object.
    constexpr auto operator()(seqan3::shape const & shape) const
    {
        return seqan3::detail::adaptor_from_functor{*this, shape, shape.size()};
    }
    //!\brief Store the shape and the window_size and return a range adaptor closure object.
    constexpr auto operator()(seqan3::shape const & shape, uint32_t const & window_size) const
    {
        return seqan3::detail::adaptor_from_functor{*this, shape, window_size};
    }
    //!\brief Store the shape and the window_size and return a range adaptor closure object.
    constexpr auto operator()(seqan3::shape const & shape, uint32_t const & window_size, uint64_t const & seed) const
    {
        return seqan3::detail::adaptor_from_functor{*this, shape, window_size, seed};
    }

    /*!\brief            Call the view's constructor with the underlying view, a seqan3::shape and a window size as
     *                   argument.
     * \param[in] urange The input range to process. Must model std::ranges::viewable_range and the reference type
     *                   of the range must model seqan3::semialphabet.
     * \param[in] shape_ The seqan3::shape to use for hashing.
     * \throws std::invalid_argument if resulting hash values would be too big for a 64 bit integer.
     * \returns          A range of converted elements.
     */

    /*!\brief            Call the view's constructor with the underlying view, a seqan3::shape and a window size as
     *                   argument.
     * \param[in] urange The input range to process. Must model std::ranges::viewable_range and the reference type
     *                   of the range must model seqan3::semialphabet.
     * \param[in] shape_ The seqan3::shape to use for hashing.
     * \throws std::invalid_argument if resulting hash values would be too big for a 64 bit integer.
     * \returns          A range of converted elements.
     */
    template <std::ranges::range urng_t>
    constexpr auto operator()(urng_t && urange, seqan3::shape const & shape, uint32_t const & window_size) const
    {
        static_assert(std::ranges::viewable_range<urng_t>,
            "The range parameter to views::minimizer cannot be a temporary of a non-view range.");
        static_assert(std::ranges::forward_range<urng_t>,
            "The range parameter to views::minimizer must model std::ranges::forward_range.");

        return std::forward<urng_t>(urange) | seqan3::views::kmer_hash(shape) | seqan3::views::minimizer(shape.size(),
                                    window_size);
    }

    template <std::ranges::range urng_t>
    constexpr auto operator()(urng_t && urange, seqan3::shape const & shape, uint32_t const & window_size,
                              uint64_t const & seed) const
    {
        static_assert(std::ranges::viewable_range<urng_t>,
            "The range parameter to views::minimizer cannot be a temporary of a non-view range.");
        static_assert(std::ranges::forward_range<urng_t>,
            "The range parameter to views::minimizer must model std::ranges::forward_range.");

        return std::forward<urng_t>(urange) | seqan3::views::kmer_hash(shape) | seqan3::views::minimizer(shape.size(),
                                    window_size, seed);
    }
};
//![adaptor_def]

} // namespace seqan3::detail

namespace seqan3::views
{

/*!\name Alphabet related views
 * \{
 */

/*!\brief                    Computes minimizers for a range with a given shape, window size and seed.
 * \tparam urng_t            The type of the range being processed. See below for requirements. [template parameter is
 *                           omitted in pipe notation]
 * \param[in] urange         The range being processed. [parameter is omitted in pipe notation]
 * \param[in] shape          The seqan3::shape that determines how to compute the hash value.
 * \param[in] window_size    The swindow size to use.
 * \param[in] seed           The seed to use. Default: 0x8F3F73B5CF1C9ADE.
 * \returns                  A range of std::size_t where each value is the hash of the resp. k-mer.
 *                           See below for the properties of the returned range.
 * \ingroup views
 *
 * \details
 *
 * \attention
 * For the alphabet size \f$\sigma\f$ of the alphabet of `urange` and the number of 1s \f$s\f$ of `shape` it must hold
 * that \f$s>\frac{64}{\log_2\sigma}\f$, i.e. hashes resulting from the shape/alphabet combination can be represented
 * in an `uint64_t`.
 *
 * ### View properties
 *
 * | Concepts and traits              | `urng_t` (underlying range type)   | `rrng_t` (returned range type)   |
 * |----------------------------------|:----------------------------------:|:--------------------------------:|
 * | std::ranges::input_range         | *required*                         | *preserved*                      |
 * | std::ranges::forward_range       | *required*                         | *preserved*                      |
 * | std::ranges::bidirectional_range |                                    | *preserved*                      |
 * | std::ranges::random_access_range |                                    | *preserved*                      |
 * | std::ranges::contiguous_range    |                                    | *lost*                           |
 * |                                  |                                    |                                  |
 * | std::ranges::viewable_range      | *required*                         | *guaranteed*                     |
 * | std::ranges::view                |                                    | *guaranteed*                     |
 * | std::ranges::sized_range         |                                    | *preserved*                      |
 * | std::ranges::common_range        |                                    | *lost*                           |
 * | std::ranges::output_range        |                                    | *lost*                           |
 * | seqan3::const_iterable_range     |                                    | *preserved*                      |
 * |                                  |                                    |                                  |
 * | std::ranges::range_reference_t   | seqan3::semialphabet               | std::size_t                      |
 *
 * See the \link views views submodule documentation \endlink for detailed descriptions of the view properties.
 *
 * ### Example
 *
 *
 *
 * \hideinitializer
 */
inline auto constexpr minimizer_hash = detail::minimizer_hash_fn{};

//!\}

} // namespace seqan3::views
