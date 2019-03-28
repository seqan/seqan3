// -----------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2019, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2019, Knut Reinert & MPI für molekulare Genetik
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at: https://github.com/seqan/seqan3/blob/master/LICENSE
// -----------------------------------------------------------------------------------------------------

/*!\file
 * \brief Contains algorithms to modify seqan3::simd::simd_type.
 * \author Marcel Ehrhardt <marcel.ehrhardt AT fu-berlin.de>
 */

#pragma once

#include <utility>

#include <seqan3/alphabet/concept.hpp>
#include <seqan3/core/metafunction/range.hpp>
#include <seqan3/core/simd/concept.hpp>
#include <seqan3/core/simd/simd_traits.hpp>
#include <seqan3/std/ranges>

namespace seqan3::detail
{

//!\brief Helper function for seqan3::simd::fill.
//!\ingroup simd
template <simd_concept simd_t, size_t... I>
constexpr simd_t fill_impl(typename simd_traits<simd_t>::scalar_type const scalar, std::index_sequence<I...>)
{
    return simd_t{((void)I, scalar)...};
}

//!\brief Helper function for seqan3::simd::iota.
//!\ingroup simd
template <simd_concept simd_t, typename scalar_t, scalar_t... I>
constexpr simd_t iota_impl(scalar_t const offset, std::integer_sequence<scalar_t, I...>)
{
    return simd_t{static_cast<scalar_t>(offset + I)...};
}

/*!\brief Transforms a batch of ranges into simd types (Structure-of-Arrays) and writes them into the output iterator.
 * \tparam simd_t        The simd type; must model seqan3::simd::simd_concept.
 * \tparam output_iter_t The output iterator type; must model std::OutputIterator with `simd_t`.
 * \tparam range_t       The input range type; must model std::ranges::ForwardRange and its reference type must model
 *                       std::ranges::ViewableRange and std::ranges::ForwardRange over seqan3::Alphabet.
 * \param out     The output iterator.
 * \param seq_rng The range over the ranges to transform.
 * \ingroup simd
 *
 * \details
 *
 * This function assumes that the batch contains simd_traits<simd_t>::length many sequences which have the same
 * size. At most simd_traits<simd_t>::length characters of the sequence are transformed or less if the ranges are
 * smaller. If less sequences are provided or the sequences have different size the behaviour is undefined.
 */
template <simd_concept simd_t,
          std::OutputIterator<simd_t> output_iter_t,
          std::ranges::InputRange range_t>
//!\cond
    requires std::ranges::ForwardRange<reference_t<range_t>> &&
             std::ranges::ViewableRange<reference_t<range_t>> &&
             seqan3::Alphabet<reference_t<reference_t<range_t>>>
//!\endcond
constexpr void transform_batch_to_soa(output_iter_t out, range_t && seq_rng)
{
    // Check if the length is identical in debug mode.
    if constexpr (std::ranges::SizedRange<range_t>)
        assert(std::ranges::size(seq_rng) == static_cast<size_t>(simd_traits<simd_t>::length));

    // stack memory for the cached iterator - sentinel pairs.
    using it_pair_t = std::pair<std::ranges::iterator_t<reference_t<range_t>>,
                                std::ranges::sentinel_t<reference_t<range_t>>>;
    std::array<it_pair_t, simd_traits<simd_t>::length> iter_cache;

    // Cache the iterator, sentinel pairs of the underlying ranges.
    for (auto && [rng, idx] : std::view::zip(seq_rng, std::view::iota(0)))
        iter_cache[idx] = {std::ranges::begin(rng), std::ranges::end(rng)};

    // Iterate over the length of the array or the maximal size of the underlying range.
    // We assume all ranges have the same length.
    simd_t simd{};
    for (size_t j = 0; (j < iter_cache.size()) || (iter_cache[j].first == iter_cache[j].second); ++j, ++out)
    {
        // Fill the simd value with the ranks of the respective alphabets.
        for (size_t i = 0; i < iter_cache.size(); ++i)
        {
            simd[i] = to_rank(*iter_cache[i].first);
            ++(iter_cache[i].first);
        }
        // Store the simd vector.
        *out = simd;
    }
}

} // namespace seqan3::detail

namespace seqan3
{

inline namespace simd
{

/*!\brief Fills a seqan3::simd::simd_type vector with a scalar value.
 * \tparam    simd_t The simd type which satisfies seqan3::simd::simd_concept.
 * \param[in] scalar The scalar value to fill the seqan3::simd::simd_type vector.
 * \ingroup simd
 *
 * \details
 *
 * \include test/snippet/core/simd/fill.cpp
 */
template <simd_concept simd_t>
constexpr simd_t fill(typename simd_traits<simd_t>::scalar_type const scalar)
{
    constexpr size_t length = simd_traits<simd_t>::length;
    return detail::fill_impl<simd_t>(scalar, std::make_index_sequence<length>{});
}

/*!\brief Fills a seqan3::simd::simd_type vector with the scalar values offset, offset+1, offset+2, ...
 * \tparam    simd_t The simd type which satisfies seqan3::simd::simd_concept.
 * \param[in] offset The scalar offset to fill the seqan3::simd::simd_type vector.
 * \ingroup simd
 *
 * \details
 *
 * \include test/snippet/core/simd/iota.cpp
 */
template <simd_concept simd_t>
constexpr simd_t iota(typename simd_traits<simd_t>::scalar_type const offset)
{
    constexpr size_t length = simd_traits<simd_t>::length;
    using scalar_type = typename simd_traits<simd_t>::scalar_type;
    return detail::iota_impl<simd_t>(offset, std::make_integer_sequence<scalar_type, length>{});
}

} // inline namespace simd

} // namespace seqan3
