// SPDX-FileCopyrightText: 2006-2025 Knut Reinert & Freie Universität Berlin
// SPDX-FileCopyrightText: 2016-2025 Knut Reinert & MPI für molekulare Genetik
// SPDX-License-Identifier: BSD-3-Clause

#pragma once

#include <thread>

#include <seqan3/core/platform.hpp>

namespace seqan3::test
{

// We don't know if the system supports hyper-threading so we use only half the threads so that the
// simd benchmark is likely to run on physical cores only.
inline uint32_t get_number_of_physical_threads()
{
    uint32_t thread_count = std::thread::hardware_concurrency();
    return (thread_count == 1) ? thread_count : thread_count >> 1;
}

} // namespace seqan3::test
