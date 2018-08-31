// ==========================================================================
//                 SeqAn - The Library for Sequence Analysis
// ==========================================================================
//
// Copyright (c) 2006-2018, Knut Reinert, FU Berlin
// Copyright (c) 2016-2018, Knut Reinert & MPI Molekulare Genetik
// All rights reserved.
//
// Redistribution and use in source and binary forms, with or without
// modification, are permitted provided that the following conditions are met:
//
//     * Redistributions of source code must retain the above copyright
//       notice, this list of conditions and the following disclaimer.
//     * Redistributions in binary form must reproduce the above copyright
//       notice, this list of conditions and the following disclaimer in the
//       documentation and/or other materials provided with the distribution.
//     * Neither the name of Knut Reinert or the FU Berlin nor the names of
//       its contributors may be used to endorse or promote products derived
//       from this software without specific prior written permission.
//
// THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
// AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
// IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
// ARE DISCLAIMED. IN NO EVENT SHALL KNUT REINERT OR THE FU BERLIN BE LIABLE
// FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL
// DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR
// SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER
// CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT
// LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY
// OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH
// DAMAGE.
//
// ==========================================================================

#include <algorithm>
#include <cstring>

#include <benchmark/benchmark.h>

static void vector_copy_benchmark(benchmark::State& state) {
    std::vector<int> x = {15, 13, 12, 10};
    for (auto _ : state)
        std::vector<int> copy{x};
}

static void memcpy_benchmark(benchmark::State& state) {
    unsigned size = state.range(0);
    char* src = new char[size];
    char* dst = new char[size];

    memset(src, '-', size);
    for (auto _ : state)
        memcpy(dst, src, size);

    int64_t bytes = int64_t(state.iterations()) * int64_t(size);
    state.counters["bytes_processed"] = bytes;
    delete[] src;
    delete[] dst;
}

static void copy_benchmark(benchmark::State& state) {
    unsigned size = state.range(0);
    char* src = new char[size];
    char* dst = new char[size];

    memset(src, '-', size);
    for (auto _ : state)
        std::copy_n(src, size, dst);

    int64_t bytes = int64_t(state.iterations()) * int64_t(size);
    state.counters["bytes_processed"] = bytes;
    delete[] src;
    delete[] dst;
}

// Register the function as a benchmark
BENCHMARK(vector_copy_benchmark);

BENCHMARK(memcpy_benchmark)->Arg(8)->Arg(64)->Arg(512);
BENCHMARK(memcpy_benchmark)->Range(4, 4<<5);
BENCHMARK(copy_benchmark)->RangeMultiplier(2)->Range(4, 4<<5);

BENCHMARK_MAIN();
