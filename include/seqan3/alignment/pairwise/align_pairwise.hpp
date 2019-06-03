// -----------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2019, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2019, Knut Reinert & MPI für molekulare Genetik
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

#include <seqan3/alignment/configuration/all.hpp>
#include <seqan3/alignment/pairwise/alignment_result.hpp>
#include <seqan3/alignment/pairwise/alignment_configurator.hpp>
#include <seqan3/alignment/pairwise/execution/all.hpp>

#include <seqan3/alphabet/gap/gapped.hpp>

#include <seqan3/core/algorithm/all.hpp>
#include <seqan3/core/metafunction/basic.hpp>
#include <seqan3/range/view/persist.hpp>

#include <seqan3/std/concepts>
#include <seqan3/std/ranges>

namespace seqan3
{

/*!\brief Computes a pairwise alignment over a sequence pair or a range of sequence pairs.
 * \ingroup pairwise_alignment
 *
 * \tparam sequence_t         The type of sequence pairs. Either a type modelling seqan3::TupleLike
 *                            with a tuple size of exactly two or std::ranges::InputRange whose value type models the
 *                            before mentioned tuple concept.
 * \tparam alignment_config_t The type of the alignment configuration; must be a seqan3::configuration.
 * \param[in] seq             A tuple with two sequences or an input range over such tuples.
 * \param[in] config          The object storing the alignment configuration.
 * \return A seqan3::alignment_range.
 *
 * \details
 *
 * This function computes the pairwise alignment of a sequence pair or a range over sequence pairs.
 * The pairs can be any class template that models the seqan3::TupleLike.
 * In case of calling the alignment with a single tuple, the value is passed as a std::view::single to the range
 * aware interface. In this case the passed class must be std::CopyConstructible and both element types of the tuple
 * must model std::ranges::ViewableRange because it is copied into the std::view::single. The same
 * is true, if the input range is a temporary non-view range.
 * Accordingly, the following example wouldn't compile:
 * ```cpp
 * std::pair p{"ACGTAGC", "AGTACGACG"};
 * auto rng = align_pairwise(p, config);
 * ```
 * You can use std::tie instead as shown in the example below.
 *
 * ### Accessing the alignment results
 *
 * For each sequence pair one or more \ref seqan3::alignment_result "seqan3::alignment_result"s can be computed.
 * The seqan3::align_pairwise function returns an seqan3::alignment_range which can be used to iterate over the
 * alignments. The alignments are then computed on-demand when iterating over the results.
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
    requires detail::AlignPairwiseSingleInput<std::remove_reference_t<sequence_t>> &&
             std::CopyConstructible<std::remove_reference_t<sequence_t>> &&
             detail::is_type_specialisation_of_v<alignment_config_t, configuration>
//!\endcond
constexpr auto align_pairwise(sequence_t && seq, alignment_config_t const & config)
{
    static_assert(std::tuple_size_v<std::remove_reference_t<sequence_t>> == 2,
                  "Alignment configuration error: Expects exactly two sequences for pairwise alignments.");

    static_assert(std::ranges::ViewableRange<std::tuple_element_t<0, std::remove_reference_t<sequence_t>>> &&
                  std::ranges::ViewableRange<std::tuple_element_t<1, std::remove_reference_t<sequence_t>>>,
                  "Alignment configuration error: The tuple elements must model std::ranges::ViewableRange.");

    return align_pairwise(std::view::single(std::forward<sequence_t>(seq)), config);
}

//!\cond
template <typename sequence_t, typename alignment_config_t>
    requires detail::AlignPairwiseRangeInputConcept<sequence_t> &&
             detail::is_type_specialisation_of_v<alignment_config_t, configuration>
constexpr auto align_pairwise(sequence_t && seq, alignment_config_t const & config)
{
    // Pipe with view::persist to allow rvalue non-view ranges.
    auto seq_view = std::forward<sequence_t>(seq) | view::persist;
    // Configure the alignment algorithm.
    auto kernel = detail::alignment_configurator::configure<decltype(seq_view)>(config);
    // Create a two-way executor for the alignment.
    detail::alignment_executor_two_way exec{std::move(seq_view), kernel};
    // Return the range over the alignments.
    return alignment_range{std::move(exec)};
}
//!\endcond

} // namespace seqan3
