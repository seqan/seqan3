// SPDX-FileCopyrightText: 2006-2025 Knut Reinert & Freie Universität Berlin
// SPDX-FileCopyrightText: 2016-2025 Knut Reinert & MPI für molekulare Genetik
// SPDX-License-Identifier: BSD-3-Clause

#include <benchmark/benchmark.h>

#include <functional>
#include <ranges>
#include <vector>

#include <seqan3/utility/simd/algorithm.hpp>
#include <seqan3/utility/simd/simd.hpp>
#include <seqan3/utility/simd/views/iota_simd.hpp>

// ============================================================================
//  iota_simd_view_benchmark
// ============================================================================

template <typename index_simd_t>
struct iota_simd_view_benchmark
{
    using index_scalar_type = typename seqan3::simd::simd_traits<index_simd_t>::scalar_type;

    seqan3::detail::iota_simd_view<index_simd_t> iota_simd_view{};

    iota_simd_view_benchmark(size_t const end_index) :
        iota_simd_view{static_cast<index_scalar_type>(0), static_cast<index_scalar_type>(end_index)}
    {}

    void operator()(index_simd_t & count) const
    {
        for (auto simd_index : iota_simd_view)
            count += simd_index;
    }
};

// ============================================================================
//  transform_iota_with_simd_fill_benchmark
// ============================================================================

template <typename index_simd_t>
struct transform_iota_with_simd_fill_benchmark
{
    using iota_transform_view_type = decltype(std::views::iota(0u, static_cast<size_t>(1))
                                              | std::views::transform(seqan3::simd::fill<index_simd_t>));

    iota_transform_view_type iota_simd_view;

    transform_iota_with_simd_fill_benchmark(size_t const end_index) :
        iota_simd_view{std::views::iota(0u, end_index) | std::views::transform(seqan3::simd::fill<index_simd_t>)}
    {}

    void operator()(index_simd_t & count) const
    {
        for (auto simd_index : iota_simd_view)
            count += simd_index;
    }
};

// ============================================================================
//  pure for loop with transform
// ============================================================================

template <typename index_simd_t>
struct for_loop_with_simd_fill_benchmark
{
    size_t end_index{};

    void operator()(index_simd_t & count) const
    {
        for (size_t index = 0; index < end_index; ++index)
            count += seqan3::simd::fill<index_simd_t>(index);
    }
};

// ============================================================================
//  pure for loop with adding vector
// ============================================================================

template <typename index_simd_t>
struct for_loop_with_simd_add_benchmark
{
    size_t end_index{};

    void operator()(index_simd_t & count) const
    {
        index_simd_t simd_index{};
        for (size_t index = 0; index < end_index; ++index, simd_index += 1)
            count += simd_index;
    }
};

// ============================================================================
//  Generic benchmark function.
// ============================================================================

template <typename index_simd_t, template <typename> typename kernel_t>
void iota_simd_benchmark(benchmark::State & state)
{
    kernel_t<index_simd_t> kernel{static_cast<size_t>(state.range(0))};
    index_simd_t count{};

    for (auto _ : state)
        std::invoke(kernel, count);

    size_t total{};
    for (size_t index = 0; index < seqan3::simd::simd_traits<index_simd_t>::length; ++index)
        total += count[index];

    state.counters["total"] = total;
}

// Baseline test using for loop and simd add.
BENCHMARK_TEMPLATE(iota_simd_benchmark, seqan3::simd::simd_type_t<uint8_t>, for_loop_with_simd_add_benchmark)
    ->Arg(std::numeric_limits<uint8_t>::max());
BENCHMARK_TEMPLATE(iota_simd_benchmark, seqan3::simd::simd_type_t<uint16_t>, for_loop_with_simd_add_benchmark)
    ->Arg(std::numeric_limits<uint16_t>::max());
BENCHMARK_TEMPLATE(iota_simd_benchmark, seqan3::simd::simd_type_t<uint32_t>, for_loop_with_simd_add_benchmark)
    ->Arg(1'000'000);
BENCHMARK_TEMPLATE(iota_simd_benchmark, seqan3::simd::simd_type_t<uint64_t>, for_loop_with_simd_add_benchmark)
    ->Arg(1'000'000);

// Baseline test using for loop and simd fill.
BENCHMARK_TEMPLATE(iota_simd_benchmark, seqan3::simd::simd_type_t<uint8_t>, for_loop_with_simd_fill_benchmark)
    ->Arg(std::numeric_limits<uint8_t>::max());
BENCHMARK_TEMPLATE(iota_simd_benchmark, seqan3::simd::simd_type_t<uint16_t>, for_loop_with_simd_fill_benchmark)
    ->Arg(std::numeric_limits<uint16_t>::max());
BENCHMARK_TEMPLATE(iota_simd_benchmark, seqan3::simd::simd_type_t<uint32_t>, for_loop_with_simd_fill_benchmark)
    ->Arg(1'000'000);
BENCHMARK_TEMPLATE(iota_simd_benchmark, seqan3::simd::simd_type_t<uint64_t>, for_loop_with_simd_fill_benchmark)
    ->Arg(1'000'000);

// Test iota view in combination with transform view that converts the scalar index to a simd vector on access.
BENCHMARK_TEMPLATE(iota_simd_benchmark, seqan3::simd::simd_type_t<uint8_t>, transform_iota_with_simd_fill_benchmark)
    ->Arg(std::numeric_limits<uint8_t>::max());
BENCHMARK_TEMPLATE(iota_simd_benchmark, seqan3::simd::simd_type_t<uint16_t>, transform_iota_with_simd_fill_benchmark)
    ->Arg(std::numeric_limits<uint16_t>::max());
BENCHMARK_TEMPLATE(iota_simd_benchmark, seqan3::simd::simd_type_t<uint32_t>, transform_iota_with_simd_fill_benchmark)
    ->Arg(1'000'000);
BENCHMARK_TEMPLATE(iota_simd_benchmark, seqan3::simd::simd_type_t<uint64_t>, transform_iota_with_simd_fill_benchmark)
    ->Arg(1'000'000);

// Test seqan3::views::iota_simd.
BENCHMARK_TEMPLATE(iota_simd_benchmark, seqan3::simd::simd_type_t<uint8_t>, iota_simd_view_benchmark)
    ->Arg(std::numeric_limits<uint8_t>::max());
BENCHMARK_TEMPLATE(iota_simd_benchmark, seqan3::simd::simd_type_t<uint16_t>, iota_simd_view_benchmark)
    ->Arg(std::numeric_limits<uint16_t>::max());
BENCHMARK_TEMPLATE(iota_simd_benchmark, seqan3::simd::simd_type_t<uint32_t>, iota_simd_view_benchmark)->Arg(1'000'000);
BENCHMARK_TEMPLATE(iota_simd_benchmark, seqan3::simd::simd_type_t<uint64_t>, iota_simd_view_benchmark)->Arg(1'000'000);

// ============================================================================
//  run
// ============================================================================

BENCHMARK_MAIN();
