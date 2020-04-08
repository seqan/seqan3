// -----------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2020, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2020, Knut Reinert & MPI für molekulare Genetik
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at: https://github.com/seqan/seqan3/blob/master/LICENSE.md
// -----------------------------------------------------------------------------------------------------

/*!\file
 * \author Mitra Darvish <mitra.darvish AT fu-berlin.de>
 * \brief Provides seqan3::views::minimiser_hash.
 */

#pragma once

#include <seqan3/range/views/kmer_hash.hpp>
#include <seqan3/range/views/minimiser.hpp>
#include <seqan3/search/kmer_index/shape.hpp>

#include <seqan3/core/debug_stream.hpp>

namespace seqan3::detail
{
// ============================================================================
//  views::minimiser_hash (adaptor instance definition)
// ============================================================================

//![adaptor_def]
//!\brief views::minimiser's range adaptor object type (non-closure).
struct minimiser_hash_fn
{
    //!\brief Store the shape and return a range adaptor closure object.
    constexpr auto operator()(seqan3::shape const & shape) const
    {
        return seqan3::detail::adaptor_from_functor{*this, shape, shape.size()};
    }
    //!\brief Store the shape and the window_size and return a range adaptor closure object.
    constexpr auto operator()(seqan3::shape const & shape, uint32_t const & window_size) const
    {
        return seqan3::detail::adaptor_from_functor{*this, shape, window_size};
    }
    //!\brief Store the shape, the window_size and the seed and return a range adaptor closure object.
    constexpr auto operator()(seqan3::shape const & shape, uint32_t const & window_size, uint64_t const & seed) const
    {
        return seqan3::detail::adaptor_from_functor{*this, shape, window_size, seed};
    }

    /*!\brief            Call the view's constructor with the underlying view, a seqan3::shape and a window size as
     *                   argument.
     * \param[in] urange The input range to process. Must model std::ranges::viewable_range and the reference type
     *                   of the range must model seqan3::semialphabet.
     * \param[in] shape_ The seqan3::shape to use for hashing.
     * \throws std::invalid_argument if the size of the shape is greater than the window_size.
     * \returns          A range of converted elements.
     */
    template <std::ranges::range urng_t>
    constexpr auto operator()(urng_t && urange, seqan3::shape const & shape, uint32_t const & window_size,
                              uint64_t const & seed = 0x8F3F73B5CF1C9ADE) const
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
                                            | std::views::transform([seed] (int i) { return i ^ seed; })
                                            | seqan3::views::minimiser(window_size - shape.size() + 1);
    }
};
//![adaptor_def]

} // namespace seqan3::detail

namespace seqan3::views
{

/*!\name Alphabet related views
 * \{
 */

/*!\brief                    Computes minimisers for a range with a given shape, window size and seed.
 * \tparam urng_t            The type of the range being processed. See below for requirements. [template parameter is
 *                           omitted in pipe notation]
 * \param[in] urange         The range being processed. [parameter is omitted in pipe notation]
 * \param[in] shape          The seqan3::shape that determines how to compute the hash value.
 * \param[in] window_size    The window size to use.
 * \param[in] seed           The seed to use. Default: 0x8F3F73B5CF1C9ADE.
 * \returns                  A range of std::integral where each value is the minimiser of the resp. window.
 *                           See below for the properties of the returned range.
 * \ingroup views
 *
 * \details
 *
 * A sequence can be presented by a small number of k-mers (minimisers). For a given shape and window size all k-mers
 * are determined in the forward strand and only the lexicographically smallest k-mer is saved for
 * one window. This process is repeated over every possible window of a sequence. If consecutive windows share a
 * minimiser, it is saved only once. It might happen that a minimiser changes only slightly when sliding
 * the window over the sequence. For instance, when a minimiser starts with a repetition of A’s, then in the next window
 * it is highly likely that the minimiser will start with a repetition of A’s as well. Because it is only one A shorter,
 * depending on how long the repetition is this might go on for multiple window shifts. Saving these only slightly
 * different minimiser makes no sense because they contain no new information about the underlying sequence plus
 * sequences with a repetition of A’s will be seen as more similar to each other than they actually are.
 * As [Marçais et al.](https://doi.org/10.1093/bioinformatics/btx235) have shown, randomizing the order of the k-mers
 * can solve this problem. Therefore, a random seed is used as a default to XOR all k-mers, thereby randomzing the
 * order. The user can change the seed to any other value he or she thinks is useful. A seed of 0 is returning the
 * lexicographical order.
 *
 * \attention
 * As f the the seqan3::views::kmer_hash the alphabet size \f$\sigma\f$ of the alphabet of `urange` and the number of
 * 1s \f$s\f$ of `shape` it must hold that \f$s>\frac{64}{\log_2\sigma}\f$, i.e. hashes resulting from the
 * shape/alphabet combination can be represented in an `uint64_t`.
 *
 * ### View properties
 *
 * | Concepts and traits              | `urng_t` (underlying range type)   | `rrng_t` (returned range type)   |
 * |----------------------------------|:----------------------------------:|:--------------------------------:|
 * | std::ranges::input_range         | *required*                         | *preserved*                      |
 * | std::ranges::forward_range       | *required*                         | *preserved*                      |
 * | std::ranges::bidirectional_range |                                    | *lost*                           |
 * | std::ranges::random_access_range |                                    | *lost*                           |
 * | std::ranges::contiguous_range    |                                    | *lost*                           |
 * |                                  |                                    |                                  |
 * | std::ranges::viewable_range      | *required*                         | *guaranteed*                     |
 * | std::ranges::view                |                                    | *guaranteed*                     |
 * | std::ranges::sized_range         |                                    | *lost*                           |
 * | std::ranges::common_range        |                                    | *lost*                           |
 * | std::ranges::output_range        |                                    | *lost*                           |
 * | seqan3::const_iterable_range     |                                    | *preserved*                      |
 * |                                  |                                    |                                  |
 * | std::ranges::range_reference_t   | seqan3::semialphabet               | std::integral                    |
 *
 * See the \link views views submodule documentation \endlink for detailed descriptions of the view properties.
 *
 * ### Example
 *
 * \include test/snippet/range/views/minimiser_hash.cpp
 *
 * \hideinitializer
 */
inline constexpr auto minimiser_hash = detail::minimiser_hash_fn{};

//!\}

} // namespace seqan3::views
