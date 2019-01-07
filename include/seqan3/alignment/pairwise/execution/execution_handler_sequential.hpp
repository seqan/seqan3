// -----------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2019, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2019, Knut Reinert & MPI für molekulare Genetik
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at: https://github.com/seqan/seqan3/blob/master/LICENSE
// -----------------------------------------------------------------------------------------------------

/*!\file
 * \brief Provides seqan3::detail::execution_handler_sequential.
 * \author Rene Rahn <rene.rahn AT fu-berlin.de>
 */

#pragma once

#include <functional>

#include <seqan3/core/platform.hpp>

namespace seqan3::detail
{

/*!\brief Handles the sequential execution of alignments.
 * \ingroup execution
 */
struct execution_handler_sequential
{
public:

    /*!\name Execution
     * \{
     */
    //!\brief Invokes the passed alignment instance in a blocking manner.
    template <typename func_result_t, typename first_batch_t, typename second_batch_t>
    void execute(std::function<func_result_t(first_batch_t const &, second_batch_t const &)> func,
                 first_batch_t const & first_batch,
                 second_batch_t const & second_batch,
                 std::function<void(decltype(func(first_batch, second_batch)))> delegate)
    {
        delegate(func(first_batch, second_batch));
    }
    //!\}
};

} // namespace seqan3
