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
#include <seqan3/std/concepts>

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
    /*!\brief Invokes the passed alignment instance in a blocking manner.
     * \tparam fn_type           The callable that needs to be invoked; must model std::Invocable with first_range_type
     *                           and second_range_type.
     * \tparam first_range_type  The type of the first range.
     * \tparam second_range_type The type of the second range.
     * \tparam delegate_type     The type of the callable invoked on the std::invoke_result of `fn_type`; must model
     *                           std::Invocable.
     *
     * \param[in] func         The callable invoking the alignment algorithm.
     * \param[in] first_range  The first range.
     * \param[in] second_range The second range.
     * \param[in] delegate     The callable invoked with the result of the alignment.
     */
    template <typename fn_type, typename first_range_type, typename second_range_type, typename delegate_type>
    //!\cond
        requires std::Invocable<fn_type, first_range_type, second_range_type> &&
                 std::Invocable<delegate_type, std::invoke_result_t<fn_type, first_range_type, second_range_type>>
    //!\endcond
    void execute(fn_type && func,
                 first_range_type && first_range,
                 second_range_type && second_range,
                 delegate_type && delegate)
    {
        delegate(func(std::forward<first_range_type>(first_range), std::forward<second_range_type>(second_range)));
    }
    //!\}
};

} // namespace seqan3
