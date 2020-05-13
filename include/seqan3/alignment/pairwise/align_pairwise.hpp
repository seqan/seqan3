// -----------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2020, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2020, Knut Reinert & MPI für molekulare Genetik
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at: https://github.com/seqan/seqan3/blob/master/LICENSE.md
// -----------------------------------------------------------------------------------------------------

/*!\file
 * \brief Provides pairwise alignment function.
 * \author Rene Rahn <rene.rahn AT fu-berlin.de>
 */

#pragma once

#include <functional>
#include <iostream>
#include <tuple>
#include <type_traits>

#include <meta/meta.hpp>

#include <seqan3/alignment/pairwise/alignment_result.hpp>
#include <seqan3/alignment/pairwise/alignment_configurator.hpp>
#include <seqan3/alignment/pairwise/detail/concept.hpp>
#include <seqan3/alignment/pairwise/detail/type_traits.hpp>
#include <seqan3/core/parallel/execution.hpp>
#include <seqan3/core/simd/simd_traits.hpp>
#include <seqan3/core/simd/simd.hpp>
#include <seqan3/core/type_traits/basic.hpp>
#include <seqan3/range/views/persist.hpp>
#include <seqan3/std/concepts>
#include <seqan3/std/ranges>

namespace seqan3
{

/*!\brief Computes the pairwise alignment for a pair of sequences or a range over sequence pairs.
 * \ingroup pairwise_alignment
 * \tparam sequence_t         The type of sequence pairs (see details for more information on the type constraints).
 * \tparam alignment_config_t The type of the alignment configuration; must be a seqan3::configuration.
 * \param[in] seq             A tuple with two sequences or a range over such tuples.
 * \param[in] config          The object storing the alignment configuration.
 * \return A seqan3::alignment_range.
 *
 * \details
 *
 * This function computes the pairwise alignment for the given sequences. During the setup phase the most efficient
 * implementation is selected depending on the configurations stored in the given seqan3::configuration object.
 * The configuration also holds settings for parallel or vectorised execution.
 *
 * ### Compute a single alignment
 *
 * In cases where only a single alignment is to be computed, the two sequences can be passed as a pair.
 * The pair can be any class template that models the seqan3::tuple_like concept. The tuple elements must model
 * std::ranges::viewable_range and std::copy_constructible.
 * Accordingly, the following example wouldn't compile:
 * ```cpp
 * std::pair p{"ACGTAGC", "AGTACGACG"};
 * auto rng = align_pairwise(p, config);
 * ```
 * You can use std::tie instead as shown in the example below.
 *
 * ### Compute multiple alignments
 *
 * In many cases one needs to compute multiple pairwise alignments. Accordingly, the align_pairwise interface allows
 * to pass a range over sequence pairs. The alignment algorithm will be configured only once for all submitted
 * alignments and then computes the alignments sequentially or in parallel depending on the given configuration.
 * Since there is always a certain amount of initial setup involving runtime checks required, it is advisable to pass
 * many sequence-pairs to this algorithm instead of repeatedly calling it with a single pair.
 *
 * ### Accessing the alignment results
 *
 * For each sequence pair one or more \ref seqan3::alignment_result "seqan3::alignment_result"s can be computed.
 * The seqan3::align_pairwise function returns an seqan3::alignment_range which can be used to iterate over the
 * alignments. If the `vectorise` configurations are omitted the alignments are computed on-demand when iterating over
 * the results. In case of a parallel execution all alignments are computed at once in parallel when calling `begin` on
 * the associated seqan3::alignment_range.
 *
 * The following snippets demonstrate the single element and the range based interface.
 *
 * \include test/snippet/alignment/pairwise/align_pairwise.cpp
 *
 * ### Exception
 *
 * Strong exception guarantee.
 *
 * Might throw std::bad_alloc if it fails to allocate the alignment matrix or seqan3::invalid_alignment_configuration
 * if the configuration is invalid.
 *
 * ### Complexity
 *
 * The complexity depends on the configured algorithm. For the \ref seqan3::align_cfg::edit "edit distance"
 * the following worst case over two input sequences of size \f$ N \f$ can be assumed:
 *
 * | Computing        | Runtime           | Space            |
 * |------------------|-------------------|------------------|
 * | score            | \f$ O(N^2/w) \f$  | \f$ O(w) \f$     |
 * | back coordinate  | \f$ O(N^2/w) \f$  | \f$ O(w) \f$     |
 * | front coordinate | \f$ O(N^2/w) \f$  | \f$ O(N^2/w) \f$ |
 * | alignment        | \f$ O(N^2/w) \f$  | \f$ O(N^2/w) \f$ |
 *
 * \f$ w \f$ is the size of a machine word.
 *
 * For all other algorithms that compute the standard dynamic programming algorithm the following worst case holds:
 *
 * | Computing        | Runtime         | Space          |
 * |------------------|-----------------|----------------|
 * | score            | \f$ O(N^2) \f$  | \f$ O(N) \f$   |
 * | back coordinate  | \f$ O(N^2) \f$  | \f$ O(N) \f$   |
 * | front coordinate | \f$ O(N^2) \f$  | \f$ O(N^2) \f$ |
 * | alignment        | \f$ O(N^2) \f$  | \f$ O(N^2) \f$ |
 *
 * In the banded case the worst case is modified as follows:
 *
 * | Computing        | Runtime         | Space          |
 * |------------------|-----------------|----------------|
 * | score            | \f$ O(N*k) \f$  | \f$ O(k) \f$   |
 * | back coordinate  | \f$ O(N*k) \f$  | \f$ O(k) \f$   |
 * | front coordinate | \f$ O(N*k) \f$  | \f$ O(N*k) \f$ |
 * | alignment        | \f$ O(N*k) \f$  | \f$ O(N*k) \f$ |
 *
 * \f$ k \f$ is the size of the band.
 *
 * ### Thread safety
 *
 * This function is re-entrant, i.e. it is always safe to call in parallel with different inputs. It is thread-safe,
 * i.e. it is safe to call in parallel with the same input under the condition that the input sequences do not change
 * when being iterated over.
 */
template <typename sequence_t, typename alignment_config_t>
//!\cond
    requires detail::align_pairwise_single_input<std::remove_reference_t<sequence_t>> &&
             std::copy_constructible<std::remove_reference_t<sequence_t>> &&
             detail::is_type_specialisation_of_v<alignment_config_t, configuration>
//!\endcond
constexpr auto align_pairwise(sequence_t && seq, alignment_config_t const & config)
{
    static_assert(std::tuple_size_v<std::remove_reference_t<sequence_t>> == 2,
                  "Alignment configuration error: Expects exactly two sequences for pairwise alignments.");

    static_assert(std::ranges::viewable_range<std::tuple_element_t<0, std::remove_reference_t<sequence_t>>> &&
                  std::ranges::viewable_range<std::tuple_element_t<1, std::remove_reference_t<sequence_t>>>,
                  "Alignment configuration error: The tuple elements must model std::ranges::viewable_range.");

    return align_pairwise(std::views::single(std::forward<sequence_t>(seq)), config);
}

//!\cond
template <typename sequence_t, typename alignment_config_t>
    requires detail::align_pairwise_range_input<sequence_t> &&
             detail::is_type_specialisation_of_v<alignment_config_t, configuration>
constexpr auto align_pairwise(sequence_t && sequences,
                              alignment_config_t const & config)
{
    using first_seq_t  = std::tuple_element_t<0, std::ranges::range_value_t<sequence_t>>;
    using second_seq_t = std::tuple_element_t<1, std::ranges::range_value_t<sequence_t>>;

    static_assert(std::ranges::random_access_range<first_seq_t> && std::ranges::sized_range<first_seq_t>,
                  "Alignment configuration error: The sequence must model random_access_range and sized_range.");
    static_assert(std::ranges::random_access_range<second_seq_t> && std::ranges::sized_range<second_seq_t>,
                  "Alignment configuration error: The sequence must model random_access_range and sized_range.");

    // Pipe with views::persist to allow rvalue non-view ranges.
    auto seq_view = std::forward<sequence_t>(sequences) | views::persist;
    // Configure the alignment algorithm.
    auto && [algorithm, complete_config] = detail::alignment_configurator::configure<decltype(seq_view)>(config);

    using traits_t = detail::alignment_configuration_traits<remove_cvref_t<decltype(complete_config)>>;
    //!brief Lambda function to translate specified parallel and or vectorised configurations into their execution rules.
    constexpr auto get_execution_rule = [] ()
    {
        if constexpr (traits_t::is_parallel)
            return seqan3::par;
        else
            return seqan3::seq;
    };

    using alignment_result_t = typename traits_t::alignment_result_type;

    auto indexed_sequence_chunk_view = views::zip(seq_view, std::views::iota(0))
                                     | views::chunk(traits_t::alignments_per_vector);

    // Create a two-way executor for the alignment.
    detail::alignment_executor_two_way executor{indexed_sequence_chunk_view,
                                                std::move(algorithm),
                                                alignment_result_t{},
                                                get_execution_rule()};
    // Return the range over the alignments.
    return alignment_range{std::move(executor)};
}
//!\endcond

} // namespace seqan3
