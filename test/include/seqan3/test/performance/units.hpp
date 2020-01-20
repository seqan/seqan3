// -----------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2020, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2020, Knut Reinert & MPI für molekulare Genetik
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at: https://github.com/seqan/seqan3/blob/master/LICENSE.md
// -----------------------------------------------------------------------------------------------------

/*!\file
 * \brief Provides predefined custom units for google benchmark.
 * \author Marcel Ehrhardt <marcel.ehrhardt AT fu-berlin.de>
 */

#pragma once

#include <benchmark/benchmark.h>

#include <seqan3/core/concept/tuple.hpp>
#include <seqan3/core/type_traits/range.hpp>
#include <seqan3/core/platform.hpp>
#include <seqan3/std/ranges>

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

//!\brief Calculates the number of cell updates for given sequences for a specific alignment config.
template <typename sequences_range_t>
inline size_t pairwise_cell_updates(sequences_range_t const & sequences_range, auto && /*align_cfg*/)
{
    size_t matrix_cells = 0u;
    for (auto && [seq1, seq2]: sequences_range)
        matrix_cells += (std::ranges::size(seq1) + 1) * (std::ranges::size(seq2) + 1);
    return matrix_cells;
}

/*!\brief This returns a counter which represents how many cell updates were done for a matrix.
 *
 * \param  cells The total number of cells processed of a complete benchmark run.
 * \return       Returns a benchmark Counter which represents CUPS (cell updates per second).
 */
inline benchmark::Counter cell_updates_per_second(size_t cells)
{
    return benchmark::Counter(cells,
                              benchmark::Counter::kIsIterationInvariantRate,
                              benchmark::Counter::OneK::kIs1000);
}

} // namespace seqan3::test
