// -----------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2019, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2019, Knut Reinert & MPI für molekulare Genetik
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at: https://github.com/seqan/seqan3/blob/master/LICENSE
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

#include <range/v3/view/single.hpp>

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
 * \tparam sequence_t         The type of sequence pairs. Either a type modelling seqan3::tuple_like_concept
 *                            with a tuple size of exactly two or std::ranges::InputRange whose value type models the
 *                            before mentioned tuple concept.
 * \tparam alignment_config_t The type of the alignment configuration; must be a seqan3::configuration.
 * \param[in] seq             A tuple with two sequences or an input range over such tuples.
 * \param[in] config          A configuration object containing \ref alignment_configuration
 *                            "alignment configurations".
 * \return Depending on the configuration, the return type of align_pairwise is either `void` if the alignment was
 * configured to call a delegate on each hit or a seqan3::alignment_range over the alignment results.
 *
 * \details
 *
 * This function computes the pairwise alignment of a sequence pair or a range of sequence pairs.
 * The type of the elements stored in the tuples must model std::ranges::ForwardRange.
 * The value category of these elements is either an lvalue or if the passed tuple or range is an lvalue, xvalues are
 * also allowed as tuple elements.
 *
 * For each alignment one or more alignment results can be computed depending on the configuration.
 * Before computing the alignment the respective alignment algorithm is configured depending on the given alignment
 * configuration. If the configuration enables a two-way execution (the default case), a seqan3::alignment_range is
 * returned. The alignment is first computed when calling begin on the returned range or the respective iterator is
 * incremented.
 *
 * ### Example
 *
 * The following snippets demonstrate the single element and the range based interface.
 *
 * \include test/snippet/alignment/pairwise/align_pairwise.cpp
 *
 * ### Exception
 *
 * Might throw std::bad_alloc if it fails to allocate the alignment matrix or seqan3::invalid_alignment_configuration
 * if the configuration is invalid.
 *
 * ### Complexity
 *
 * The complexity depends on the configured algorithm. For the \ref seqan3::align_cfg::edit "edit distance"
 * the following worst cases over two input sequences of size \f$ N \f$ can be assumed:
 *
 * | Component        | Runtime           | Space            |
 * |------------------|-------------------|------------------|
 * | Score only       | \f$ O(N*N/w) \f$  | \f$ O(w) \f$     |
 * | Back coordinate  | \f$ O(N*N/w) \f$  | \f$ O(w) \f$     |
 * | Begin coordinate | \f$ O(N*N/w) \f$  | \f$ O(N*N/w) \f$ |
 * | Alignment        | \f$ O(N*N/w) \f$  | \f$ O(N*N/w) \f$ |
 *
 * \f$ w \f$ is the size of a machine word.
 *
 * For all other algorithms that compute the standard dynamic programming algorithm the following worst cases hold:
 *
 * | Component        | Runtime         | Space          |
 * |------------------|-----------------|----------------|
 * | Score only       | \f$ O(N^2) \f$  | \f$ O(N) \f$   |
 * | Back coordinate  | \f$ O(N^2) \f$  | \f$ O(N) \f$   |
 * | Begin coordinate | \f$ O(N^2) \f$  | \f$ O(N^2) \f$ |
 * | Alignment        | \f$ O(N^2) \f$  | \f$ O(N^2) \f$ |
 *
 * ### Thread safety
 *
 * As long as no shared state is passed to the algorithm via the configuration object, the algorithm is thread safe and
 * can be invoked in a concurrent environment.
 */
template <typename sequence_t,
          typename alignment_config_t>
//!\cond
    requires tuple_like_concept<sequence_t> &&
             detail::is_type_specialisation_of_v<remove_cvref_t<alignment_config_t>, configuration>
//!\endcond
constexpr auto align_pairwise(sequence_t && seq, alignment_config_t && config)
{
    static_assert(std::tuple_size_v<std::remove_reference_t<sequence_t>> == 2,
                  "Expects exactly two sequences for pairwise alignments.");

    return align_pairwise(ranges::view::single(std::forward<sequence_t>(seq)) | std::view::common,
                          std::forward<alignment_config_t>(config));
}

//!\copydoc align_pairwise
template <std::ranges::InputRange sequence_t, typename alignment_config_t>
//!\cond
    requires detail::is_type_specialisation_of_v<remove_cvref_t<alignment_config_t>, configuration> &&
             tuple_like_concept<value_type_t<std::ranges::iterator_t<std::remove_reference_t<sequence_t>>>>
//!\endcond
constexpr auto align_pairwise(sequence_t && seq, alignment_config_t const & config)
{
    // Wrap in persist to make the code also work with temporaries that are not views.
    auto seq_view = view::persist(std::forward<sequence_t>(seq));
    // Configure the alignment algorithm.
    auto kernel = detail::alignment_configurator::configure(seq_view, config);
    // Create a two-way executor for the alignment.
    detail::alignment_executor_two_way exec{std::move(seq_view), kernel};
    // Return the range over the alignments.
    return alignment_range{std::move(exec)};
}

} // namespace seqan3
