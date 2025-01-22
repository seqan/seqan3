// SPDX-FileCopyrightText: 2006-2025 Knut Reinert & Freie Universität Berlin
// SPDX-FileCopyrightText: 2016-2025 Knut Reinert & MPI für molekulare Genetik
// SPDX-License-Identifier: BSD-3-Clause

#include <benchmark/benchmark.h>

#include <array>
#include <cstdlib>

#include <seqan3/utility/simd/algorithm.hpp>
#include <seqan3/utility/simd/simd.hpp>
#include <seqan3/utility/simd/simd_traits.hpp>

// ----------------------------------------------------------------------------
// Helper functions
// ----------------------------------------------------------------------------

inline auto make_matrix()
{
    using simd_t = seqan3::simd::simd_type_t<int8_t>;

    std::array<simd_t, seqan3::simd::simd_traits<simd_t>::length> matrix;
    for (size_t i = 0; i < matrix.size(); ++i)
        for (size_t j = 0; j < matrix.size(); ++j)
            matrix[i][j] = std::rand() % 10;

    return matrix;
}

template <typename simd_t>
inline auto reduce(simd_t const & vec)
{
    size_t sum = 0;
    for (size_t i = 0; i < seqan3::simd::simd_traits<simd_t>::length; ++i)
        sum += vec[i];

    return sum;
}

// ----------------------------------------------------------------------------
// Benchhmark transpose
// ----------------------------------------------------------------------------

static void transpose(benchmark::State & state)
{
    size_t sum = 0;

    auto matrix = make_matrix();

    for (auto _ : state)
    {
        for (size_t i = 0; i < 100; ++i)
        {
            seqan3::simd::transpose(matrix);

            state.PauseTiming();
            sum += reduce(matrix[std::rand() % matrix.size()]);
            state.ResumeTiming();
        }
    }

    state.counters["checksum"] = sum;
}

BENCHMARK(transpose);

template <typename source_t, typename target_t>
static void upcast(benchmark::State & state)
{
    source_t src = seqan3::simd::iota<source_t>(std::rand() % 100);
    target_t target{};
    size_t sum = 0;

    for (auto _ : state)
    {
        for (size_t i = 0; i < 1000; ++i)
        {
            target = seqan3::simd::upcast<target_t>(src);

            state.PauseTiming();
            sum += reduce(target);
            state.ResumeTiming();
        }
    }

    state.counters["checksum"] = sum;
}

// ----------------------------------------------------------------------------
// Benchhmark upcast
// ----------------------------------------------------------------------------

BENCHMARK_TEMPLATE(upcast, seqan3::simd::simd_type_t<int8_t>, seqan3::simd::simd_type_t<int16_t>);
BENCHMARK_TEMPLATE(upcast, seqan3::simd::simd_type_t<int8_t>, seqan3::simd::simd_type_t<int32_t>);
BENCHMARK_TEMPLATE(upcast, seqan3::simd::simd_type_t<int8_t>, seqan3::simd::simd_type_t<int64_t>);
BENCHMARK_TEMPLATE(upcast, seqan3::simd::simd_type_t<int16_t>, seqan3::simd::simd_type_t<int32_t>);
BENCHMARK_TEMPLATE(upcast, seqan3::simd::simd_type_t<int16_t>, seqan3::simd::simd_type_t<int64_t>);
BENCHMARK_TEMPLATE(upcast, seqan3::simd::simd_type_t<int32_t>, seqan3::simd::simd_type_t<int64_t>);

BENCHMARK_MAIN();
