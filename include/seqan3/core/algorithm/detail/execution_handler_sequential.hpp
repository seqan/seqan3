// -----------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2020, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2020, Knut Reinert & MPI für molekulare Genetik
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at: https://github.com/seqan/seqan3/blob/master/LICENSE.md
// -----------------------------------------------------------------------------------------------------

/*!\file
 * \brief Provides seqan3::detail::execution_handler_sequential.
 * \author Rene Rahn <rene.rahn AT fu-berlin.de>
 */

#pragma once

#include <functional>

#include <seqan3/alignment/pairwise/detail/concept.hpp>
#include <seqan3/core/platform.hpp>
#include <seqan3/range/views/type_reduce.hpp>
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

    /*!\brief Takes underlying range of sequence pairs and invokes an alignment on each instance.
     * \tparam algorithm_t              The type of the alignment algorithm.
     * \tparam indexed_sequence_pairs_t The type of underlying sequence pairs annotated with an index.
     * \tparam delegate_type            The type of the callable invoked on the std::invoke_result of `algorithm_t`.
     *
     * \param[in] algorithm              The alignment algorithm to invoke.
     * \param[in] indexed_sequence_pairs The range of underlying annotated sequence pairs to be aligned.
     * \param[in] delegate               A callable which will be invoked on each result of the computed alignments.
     */
    template <typename algorithm_t, typename indexed_sequence_pairs_t, typename delegate_type>
    void execute(algorithm_t && algorithm,
                 indexed_sequence_pairs_t && indexed_sequence_pairs,
                 delegate_type && delegate)
    {
        algorithm(std::forward<indexed_sequence_pairs_t>(indexed_sequence_pairs),
                  std::forward<delegate_type>(delegate));
    }

    //!\brief Waits for the submitted alignments jobs to finish. (Noop).
    void wait() noexcept
    {
        // noop
    }
};

} // namespace seqan3
