// -----------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2019, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2019, Knut Reinert & MPI für molekulare Genetik
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at: https://github.com/seqan/seqan3/blob/master/LICENSE
// -----------------------------------------------------------------------------------------------------

/*!\file
 * \brief Contains predefined custom units for google benchmark.
 * \author Marcel Ehrhardt <marcel.ehrhardt AT fu-berlin.de>
 */

#pragma once

#include <benchmark/benchmark.h>

#include <seqan3/core/platform.hpp>

namespace seqan3::test
{

/*!\brief This returns a counter which represents how many bytes were processed per second.
 *
 * \param  bytes The total number of bytes processed of a complete benchmark run.
 * \return       Returns a benchmark Counter which represents bytes/s.
 */
inline benchmark::Counter bytes_per_second(size_t bytes)
{
    return benchmark::Counter(bytes,
                              benchmark::Counter::kIsIterationInvariantRate,
                              benchmark::Counter::OneK::kIs1024);
}

} // namespace seqan3::test
