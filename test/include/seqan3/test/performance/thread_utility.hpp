// -----------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2020, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2020, Knut Reinert & MPI für molekulare Genetik
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at: https://github.com/seqan/seqan3/blob/master/LICENSE.md
// -----------------------------------------------------------------------------------------------------

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
