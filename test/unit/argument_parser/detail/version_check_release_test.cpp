// -----------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2019, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2019, Knut Reinert & MPI für molekulare Genetik
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at: https://github.com/seqan/seqan3/blob/master/LICENSE
// -----------------------------------------------------------------------------------------------------

#if !defined(NDEBUG)
    #define NDEBUG // test in release mode
#endif

#include <string>

// Define global app name component to disambiguate from release test
// when tests are run in parallel.
static const std::string app_name_component{"release"};

#include "version_check_test.hpp"
