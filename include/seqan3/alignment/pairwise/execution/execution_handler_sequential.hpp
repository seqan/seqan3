// -----------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2019, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2019, Knut Reinert & MPI für molekulare Genetik
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at: https://github.com/seqan/seqan3/blob/master/LICENSE.md
// -----------------------------------------------------------------------------------------------------

/*!\file
 * \brief Provides seqan3::detail::execution_handler_sequential.
 * \author Rene Rahn <rene.rahn AT fu-berlin.de>
 */

#pragma once

#include <functional>

#include <seqan3/core/platform.hpp>
#include <seqan3/range/views/view_all.hpp>
#include <seqan3/std/concepts>
#include <seqan3/std/ranges>

namespace seqan3::detail
{

/*!\brief Handles the sequential execution of alignments.
 * \ingroup execution
 */
struct execution_handler_sequential
{
public:

    /*!\brief Invokes the passed alignment instance in a blocking manner.
     * \tparam algorithm_t       The callable that needs to be invoked; must model std::invocable with first_range_type
     *                           and second_range_type.
     * \tparam first_range_type  The type of the first range; must model std::ranges::view.
     * \tparam second_range_type The type of the second range; must model std::ranges::view.
     * \tparam delegate_type     The type of the callable invoked on the std::invoke_result of `algorithm_t`; must model
     *                           std::invocable.
     *
     * \param[in] algorithm    The alignment algorithm to invoke.
     * \param[in] idx          The index of the current processed sequence pair.
     * \param[in] first_range  The first range.
     * \param[in] second_range The second range.
     * \param[in] delegate     The callable invoked with the result of the alignment.
     */
    template <typename algorithm_t, typename first_range_type, typename second_range_type, typename delegate_type>
    //!\cond
        requires std::invocable<algorithm_t, size_t const, first_range_type, second_range_type> &&
                 std::invocable<delegate_type, std::invoke_result_t<algorithm_t,
                                                                    size_t const,
                                                                    first_range_type,
                                                                    second_range_type>>
    //!\endcond
    void execute(algorithm_t && algorithm,
                 size_t const idx,
                 first_range_type first_range,
                 second_range_type second_range,
                 delegate_type && delegate)
    {
        static_assert(std::ranges::view<first_range_type>, "Expected a view!");
        static_assert(std::ranges::view<second_range_type>, "Expected a view!");

        delegate(algorithm(idx, std::move(first_range), std::move(second_range)));
    }

    /*!\brief Takes underlying range of sequence pairs and invokes an alignment on each instance.
     * \tparam algorithm_t              The type of the alignment algorithm.
     * \tparam indexed_sequence_pairs_t The type of underlying sequence pairs annotated with an index;
     *                                  must model std::ranges::forward_range.
     * \tparam delegate_type            The type of the callable invoked on the std::invoke_result of `algorithm_t`.
     *
     * \param[in] algorithm              The alignment algorithm to invoke.
     * \param[in] indexed_sequence_pairs The range of underlying annotated sequence pairs to be aligned.
     * \param[in] delegate               A callable which will be invoked on each result of the computed alignments.
     */
    template <typename algorithm_t, std::ranges::forward_range indexed_sequence_pairs_t, typename delegate_type>
    void execute(algorithm_t && algorithm,
                 indexed_sequence_pairs_t indexed_sequence_pairs,
                 delegate_type && delegate)
    {
        using std::get;
        for (auto && [sequence_pair, idx] : indexed_sequence_pairs)
            execute(std::forward<algorithm_t>(algorithm), idx, get<0>(sequence_pair), get<1>(sequence_pair), delegate);
    }

    //!\brief Waits for the submitted alignments jobs to finish. (Noop).
    void wait() noexcept
    {
        // noop
    }
};

} // namespace seqan3
