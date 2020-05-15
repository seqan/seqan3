// -----------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2020, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2020, Knut Reinert & MPI für molekulare Genetik
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at: https://github.com/seqan/seqan3/blob/master/LICENSE.md
// -----------------------------------------------------------------------------------------------------

#include <seqan3/std/ranges>
#include <vector>

#include <benchmark/benchmark.h>

#include <seqan3/core/simd/simd_algorithm.hpp>
#include <seqan3/core/simd/simd.hpp>
#include <seqan3/core/simd/view_iota_simd.hpp>

// ============================================================================
//  simd_iota_iterator
// ============================================================================

template <typename simd_t>
struct simd_iota_iterator_function
{
    using scalar_type = typename seqan3::simd::simd_traits<simd_t>::scalar_type;

    seqan3::views::iota_simd<simd_t> simd_iota_view;

    simd_iota_iterator_function(size_t const end_index) :
        simd_iota_view{static_cast<scalar_type>(0), static_cast<scalar_type>(end_index)}
    {}

    void operator()(simd_t & count) const
    {
        for (auto simd_index : simd_iota_view)
            count += simd_index;
    }
};

// ============================================================================
//  iota_transform
// ============================================================================

template <typename simd_t>
struct transform_to_simd
{
    template <typename index_t>
    simd_t operator()(index_t const index) const
    {
        return seqan3::simd::fill<simd_t>(index);
    }
};

template <typename simd_t>
struct iota_transform_function
{
    using iota_transform_view_type = decltype(std::views::iota(0u, static_cast<size_t>(1))
                                            | std::views::transform(transform_to_simd<simd_t>{}));

    iota_transform_view_type iota_simd_view;

    iota_transform_function(size_t const end_index) :
        iota_simd_view{std::views::iota(0u, end_index) | std::views::transform(transform_to_simd<simd_t>{})}
    {}

    void operator()(simd_t & count) const
    {
        for (auto simd_index : iota_simd_view)
            count += simd_index;
    }
};

// ============================================================================
//  pure for loop with transform
// ============================================================================

template <typename simd_t>
struct for_loop_with_simd_fill
{
    size_t end_index{};

    void operator()(simd_t & count) const
    {
        for (size_t index = 0; index < end_index; ++index)
            count += transform_to_simd<simd_t>{}(index);
    }
};

// ============================================================================
//  pure for loop with adding vector
// ============================================================================

template <typename simd_t>
struct for_loop_with_simd_add
{
    size_t end_index{};

    void operator()(simd_t & count) const
    {
        simd_t simd_index{};
        for (size_t index = 0; index < end_index; ++index, ++simd_index)
            count += simd_index;
    }
};

// ============================================================================
//  Generic benchmark function.
// ============================================================================

template <typename simd_t, template <typename> typename kernel_t>
void simd_iota_benchmark(benchmark::State& state)
{
    kernel_t<simd_t> kernel{static_cast<size_t>(state.range(0))};
    simd_t count{};

    for (auto _ : state)
        std::invoke(kernel, count);

    size_t total{};
    for (size_t index = 0; index < seqan3::simd::simd_traits<simd_t>::length; ++index)
        total += count[index];

    state.counters["total"] = total;
}

// Baseline test using for loop and simd add.
BENCHMARK_TEMPLATE(simd_iota_benchmark,
                   seqan3::simd::simd_type_t<uint8_t>,
                   for_loop_with_simd_add)->Arg(std::numeric_limits<uint8_t>::max());
BENCHMARK_TEMPLATE(simd_iota_benchmark,
                   seqan3::simd::simd_type_t<uint16_t>,
                   for_loop_with_simd_add)->Arg(std::numeric_limits<uint16_t>::max());
BENCHMARK_TEMPLATE(simd_iota_benchmark,
                   seqan3::simd::simd_type_t<uint32_t>,
                   for_loop_with_simd_add)->Arg(1'000'000);
BENCHMARK_TEMPLATE(simd_iota_benchmark,
                   seqan3::simd::simd_type_t<uint64_t>,
                   for_loop_with_simd_add)->Arg(1'000'000);

// Baseline test using for loop and simd fill.
BENCHMARK_TEMPLATE(simd_iota_benchmark,
                   seqan3::simd::simd_type_t<uint8_t>,
                   for_loop_with_simd_fill)->Arg(std::numeric_limits<uint8_t>::max());
BENCHMARK_TEMPLATE(simd_iota_benchmark,
                   seqan3::simd::simd_type_t<uint16_t>,
                   for_loop_with_simd_fill)->Arg(std::numeric_limits<uint16_t>::max());
BENCHMARK_TEMPLATE(simd_iota_benchmark,
                   seqan3::simd::simd_type_t<uint32_t>,
                   for_loop_with_simd_fill)->Arg(1'000'000);
BENCHMARK_TEMPLATE(simd_iota_benchmark,
                   seqan3::simd::simd_type_t<uint64_t>,
                   for_loop_with_simd_fill)->Arg(1'000'000);

// Test iota view in combination with transform view that convers the scalar index to a simd vector on access.
BENCHMARK_TEMPLATE(simd_iota_benchmark,
                   seqan3::simd::simd_type_t<uint8_t>,
                   iota_transform_function)->Arg(std::numeric_limits<uint8_t>::max());
BENCHMARK_TEMPLATE(simd_iota_benchmark,
                   seqan3::simd::simd_type_t<uint16_t>,
                   iota_transform_function)->Arg(std::numeric_limits<uint16_t>::max());
BENCHMARK_TEMPLATE(simd_iota_benchmark,
                   seqan3::simd::simd_type_t<uint32_t>,
                   iota_transform_function)->Arg(1'000'000);
BENCHMARK_TEMPLATE(simd_iota_benchmark,
                   seqan3::simd::simd_type_t<uint64_t>,
                   iota_transform_function)->Arg(1'000'000);

// Test explicit simd iota iterator view.
BENCHMARK_TEMPLATE(simd_iota_benchmark,
                   seqan3::simd::simd_type_t<uint8_t>,
                   simd_iota_iterator_function)->Arg(std::numeric_limits<uint8_t>::max());
BENCHMARK_TEMPLATE(simd_iota_benchmark,
                   seqan3::simd::simd_type_t<uint16_t>,
                   simd_iota_iterator_function)->Arg(std::numeric_limits<uint16_t>::max());
BENCHMARK_TEMPLATE(simd_iota_benchmark,
                   seqan3::simd::simd_type_t<uint32_t>,
                   simd_iota_iterator_function)->Arg(1'000'000);
BENCHMARK_TEMPLATE(simd_iota_benchmark,
                   seqan3::simd::simd_type_t<uint64_t>,
                   simd_iota_iterator_function)->Arg(1'000'000);

// ============================================================================
//  run
// ============================================================================

BENCHMARK_MAIN();
