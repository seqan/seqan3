// SPDX-FileCopyrightText: 2006-2025 Knut Reinert & Freie Universität Berlin
// SPDX-FileCopyrightText: 2016-2025 Knut Reinert & MPI für molekulare Genetik
// SPDX-License-Identifier: BSD-3-Clause

/*!\file
 * \brief Provides pairwise alignment function.
 * \author Rene Rahn <rene.rahn AT fu-berlin.de>
 */

#pragma once

#include <concepts>
#include <functional>
#include <iostream>
#include <ranges>
#include <tuple>
#include <type_traits>

#include <seqan3/alignment/pairwise/alignment_configurator.hpp>
#include <seqan3/alignment/pairwise/alignment_result.hpp>
#include <seqan3/alignment/pairwise/detail/concept.hpp>
#include <seqan3/alignment/pairwise/detail/type_traits.hpp>
#include <seqan3/core/algorithm/algorithm_result_generator_range.hpp>
#include <seqan3/core/algorithm/detail/algorithm_executor_blocking.hpp>
#include <seqan3/utility/simd/simd.hpp>
#include <seqan3/utility/simd/simd_traits.hpp>
#include <seqan3/utility/type_traits/basic.hpp>

namespace seqan3
{

/*!\brief Computes the pairwise alignment for a pair of sequences or a range over sequence pairs.
 * \ingroup alignment_pairwise
 * \tparam sequence_t         The type of sequence pairs (see details for more information on the type constraints).
 * \tparam alignment_config_t The type of the alignment configuration; must be a seqan3::configuration.
 * \param[in] seq             A tuple with two sequences or a range over such tuples.
 * \param[in] config          The object storing the alignment configuration.
 * \return A seqan3::algorithm_result_generator_range.
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
 * The following example demonstrates how an alignment is computed from a std::pair of sequences.
 * \snippet test/snippet/alignment/pairwise/align_pairwise.cpp example1
 * Alternatively, you can use std::tie as shown in the example below:
 * \snippet test/snippet/alignment/pairwise/align_pairwise.cpp example2
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
 * The seqan3::align_pairwise function returns an seqan3::algorithm_result_generator_range which can be used to iterate
 * over the alignments. If the `vectorised` configurations are omitted the alignments are computed on-demand when
 * iterating over the results. In case of a parallel execution all alignments are computed at once in parallel when
 * calling `begin` on the associated seqan3::algorithm_result_generator_range.
 *
 * The following snippets demonstrate the single element and the range based interface.
 *
 * \snippet test/snippet/alignment/pairwise/align_pairwise.cpp start
 * \snippet test/snippet/alignment/pairwise/align_pairwise.cpp example3
 * \snippet test/snippet/alignment/pairwise/align_pairwise.cpp end
 *
 * ### Exception
 *
 * Strong exception guarantee.
 *
 * Might throw std::bad_alloc if it fails to allocate the alignment matrix or seqan3::invalid_alignment_configuration
 * if the configuration is invalid.
 * Throws std::runtime_error if seqan3::align_cfg::parallel has been specified without a `thread_count` value.
 *
 * ### Complexity
 *
 * The complexity depends on the configured algorithm. For the \ref seqan3::align_cfg::edit_scheme "edit distance"
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
    requires detail::align_pairwise_single_input<sequence_t>
          && std::copy_constructible<std::remove_reference_t<sequence_t>>
          && detail::is_type_specialisation_of_v<alignment_config_t, configuration>
constexpr auto align_pairwise(sequence_t && seq, alignment_config_t const & config)
{
    using std::get;

    if constexpr (std::is_lvalue_reference_v<sequence_t>) // Forward tuple elements as references.
    {
        return align_pairwise(std::tie(get<0>(seq), get<1>(seq)), config);
    }
    else
    {
        static_assert(std::tuple_size_v<std::remove_reference_t<sequence_t>> == 2,
                      "Alignment configuration error: Expects exactly two sequences for pairwise alignments.");

        static_assert(std::ranges::viewable_range<std::tuple_element_t<0, std::remove_reference_t<sequence_t>>>
                          && std::ranges::viewable_range<std::tuple_element_t<1, std::remove_reference_t<sequence_t>>>,
                      "Alignment configuration error: The tuple elements must model std::ranges::viewable_range.");

        return align_pairwise(std::views::single(std::forward<sequence_t>(seq)), config);
    }
}

//!\cond
template <typename sequence_t, typename alignment_config_t>
    requires detail::align_pairwise_range_input<sequence_t>
          && detail::is_type_specialisation_of_v<alignment_config_t, configuration>
constexpr auto align_pairwise(sequence_t && sequences, alignment_config_t const & config)
{
    using first_seq_t = std::tuple_element_t<0, std::ranges::range_value_t<sequence_t>>;
    using second_seq_t = std::tuple_element_t<1, std::ranges::range_value_t<sequence_t>>;

    static_assert(std::ranges::random_access_range<first_seq_t> && std::ranges::sized_range<first_seq_t>,
                  "Alignment configuration error: The sequence must model random_access_range and sized_range.");
    static_assert(std::ranges::random_access_range<second_seq_t> && std::ranges::sized_range<second_seq_t>,
                  "Alignment configuration error: The sequence must model random_access_range and sized_range.");

    // Pipe with std::views::all to allow rvalue non-view ranges.
    auto seq_view = std::forward<sequence_t>(sequences) | std::views::all;
    // Configure the alignment algorithm.
    auto && [algorithm, complete_config] = detail::alignment_configurator::configure<decltype(seq_view)>(config);

    using complete_config_t = std::remove_cvref_t<decltype(complete_config)>;
    using traits_t = detail::alignment_configuration_traits<complete_config_t>;

    auto indexed_sequence_chunk_view =
        views::zip(seq_view, std::views::iota(0)) | views::chunk(traits_t::alignments_per_vector);

    using indexed_sequences_t = decltype(indexed_sequence_chunk_view);
    using alignment_result_t = typename traits_t::alignment_result_type;
    using execution_handler_t = std::conditional_t<complete_config_t::template exists<align_cfg::parallel>(),
                                                   detail::execution_handler_parallel,
                                                   detail::execution_handler_sequential>;
    using executor_t = detail::
        algorithm_executor_blocking<indexed_sequences_t, decltype(algorithm), alignment_result_t, execution_handler_t>;

    // Select the execution handler for the alignment configuration.
    auto select_execution_handler = [parallel = complete_config.get_or(align_cfg::parallel{})]()
    {
        if constexpr (std::same_as<execution_handler_t, detail::execution_handler_parallel>)
        {
            auto thread_count = parallel.thread_count;
            if (!thread_count)
                throw std::runtime_error{"You must configure the number of threads in seqan3::align_cfg::parallel."};

            return execution_handler_t{*thread_count};
        }
        else
        {
            return execution_handler_t{};
        }
    };

    if constexpr (traits_t::is_one_way_execution) // Just compute alignment and wait until all alignments are computed.
        select_execution_handler().bulk_execute(algorithm,
                                                indexed_sequence_chunk_view,
                                                get<align_cfg::on_result>(complete_config).callback);
    else // Require two way execution: return the range over the alignments.
        return algorithm_result_generator_range{executor_t{std::move(indexed_sequence_chunk_view),
                                                           std::move(algorithm),
                                                           alignment_result_t{},
                                                           select_execution_handler()}};
}
//!\endcond

} // namespace seqan3
