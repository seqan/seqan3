// SPDX-FileCopyrightText: 2006-2025 Knut Reinert & Freie Universität Berlin
// SPDX-FileCopyrightText: 2016-2025 Knut Reinert & MPI für molekulare Genetik
// SPDX-License-Identifier: BSD-3-Clause

/*!\file
 * \brief Provides API compatibility for newer google/benchmark versions.
 * \author Enrico Seiler <enrico.seiler AT fu-berlin.de>
 */

#pragma once

#include <benchmark/benchmark.h>

#include <seqan3/core/platform.hpp>

/*!\brief Provide benchmark::Benchmark for old google/benchmark versions.
 * \sa https://github.com/google/benchmark/pull/2101
 * \details
 *
 * Google benchmark moved its `Benchmark` class from `benchmark::internal` to `benchmark`.
 * In old versions, `benchmark::Benchmark` does not exist.
 * In new versions, `benchmark::internal::Benchmark` will trigger a deprecation warning.
 *
 * Here is a reduced example of the change:
 * ```cpp
 * #ifdef OLD
 * namespace benchmark
 * {
 *     namespace internal
 *     {
 *         class Benchmark{};
 *     }
 * }
 * #endif
 *
 * #ifdef NEW
 * namespace benchmark
 * {
 *     class Benchmark{};
 *
 *     namespace internal
 *     {
 *         using Benchmark [[deprecated("Use benchmark::Benchmark instead")]] = ::benchmark::Benchmark;
 *     }
 * }
 * #endif
 * ```
 *
 * Because we do not necessarily require a specific version of benchmark, and, at the time of writing, the new change
 * is not yet included in a release, we want to already use the new API but also not break the old API.
 *
 * ### Version 1
 * https://godbolt.org/z/zzWdojW4G
 *
 * We can simply alias `benchmark::internal::Benchmark` to `benchmark::Benchmark` but ignore the deprecation warning.
 * The advantage of this version is that it is easy to understand (both code and intention), but it will break once
 * `benchmark::internal::Benchmark` is removed.
 *
 * ```cpp
 * namespace benchmark
 * {
 *
 * #pragma GCC diagnostic push
 * #pragma GCC diagnostic ignored "-Wdeprecated-declarations"
 * using Benchmark = ::benchmark::internal::Benchmark;
 * #pragma GCC diagnostic pop
 *
 * }
 * ```
 *
 * ### Version 2
 * https://godbolt.org/z/1M9f8fojv
 *
 * We use an anonymous namespace within `benchmark` to resolve `Benchmark` from either location without triggering
 * deprecation warnings. The anonymous namespace brings both `benchmark` and `benchmark::internal` into scope, then
 * aliases the first available `Benchmark` type. This works for both old and new versions and is future-proof, albeit
 * harder to grasp (both code and intention).
 *
 * ```cpp
 * namespace benchmark
 * {
 *
 * namespace
 * {
 *
 * using namespace benchmark;
 * using namespace benchmark::internal;
 * using Benchmark = Benchmark;
 *
 * }
 *
 * }
 * ```
 */
namespace benchmark
{

namespace
{

using namespace benchmark;
using namespace benchmark::internal;
using Benchmark = Benchmark;

} // namespace

} // namespace benchmark
