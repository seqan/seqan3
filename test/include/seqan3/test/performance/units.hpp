// SPDX-FileCopyrightText: 2006-2025 Knut Reinert & Freie Universität Berlin
// SPDX-FileCopyrightText: 2016-2025 Knut Reinert & MPI für molekulare Genetik
// SPDX-License-Identifier: BSD-3-Clause

/*!\file
 * \brief Provides predefined custom units for google benchmark.
 * \author Marcel Ehrhardt <marcel.ehrhardt AT fu-berlin.de>
 */

#pragma once

#include <benchmark/benchmark.h>

#include <ranges>

#include <seqan3/alignment/configuration/align_config_band.hpp>
#include <seqan3/core/platform.hpp>
#include <seqan3/core/range/type_traits.hpp>
#include <seqan3/utility/tuple/concept.hpp>

namespace seqan3::test
{

/*!\brief This returns a counter which represents how many bytes were processed per second.
 *
 * \param  bytes The total number of bytes processed of a complete benchmark run.
 * \return       Returns a benchmark Counter which represents bytes/s.
 */
inline benchmark::Counter bytes_per_second(size_t bytes)
{
    return benchmark::Counter(bytes, benchmark::Counter::kIsIterationInvariantRate, benchmark::Counter::OneK::kIs1024);
}

//!\brief Calculates the number of cell updates for given sequences for a specific alignment config.
template <typename sequences_range_t>
inline size_t pairwise_cell_updates(sequences_range_t const & sequences_range, [[maybe_unused]] auto && align_cfg)
{
    using config_t = std::remove_cvref_t<decltype(align_cfg)>;

    auto count_cells = [&](auto && seq1, auto && seq2)
    {
        size_t const columns = std::ranges::size(seq1) + 1;
        size_t const rows = std::ranges::size(seq2) + 1;

        if constexpr (config_t::template exists<seqan3::align_cfg::band_fixed_size>())
        {
            using std::get;
            auto const band_cfg = get<seqan3::align_cfg::band_fixed_size>(align_cfg);

            int32_t const lower_diagonal = band_cfg.lower_diagonal;
            int32_t const upper_diagonal = band_cfg.upper_diagonal;
            size_t matrix_cells = 0;
            for (int32_t column_id = 0; column_id < static_cast<int32_t>(columns); ++column_id)
            {
                // the position of the upper-row (inclusive)
                int32_t upper_row_id = std::clamp<int32_t>(column_id - upper_diagonal, 0, rows);
                // the position of the lower-row (exclusive)
                int32_t lower_row_id = std::clamp<int32_t>(column_id - lower_diagonal + 1, 0, rows);
                matrix_cells += lower_row_id - upper_row_id;
            }

            return matrix_cells;
        }
        else
        {
            return columns * rows;
        }
    };

    size_t matrix_cells = 0u;
    for (auto && [seq1, seq2] : sequences_range)
        matrix_cells += count_cells(seq1, seq2);

    return matrix_cells;
}

/*!\brief This returns a counter which represents how many cell updates were done for a matrix.
 *
 * \param  cells The total number of cells processed of a complete benchmark run.
 * \return       Returns a benchmark Counter which represents CUPS (cell updates per second).
 */
inline benchmark::Counter cell_updates_per_second(size_t cells)
{
    return benchmark::Counter(cells, benchmark::Counter::kIsIterationInvariantRate, benchmark::Counter::OneK::kIs1000);
}

} // namespace seqan3::test
