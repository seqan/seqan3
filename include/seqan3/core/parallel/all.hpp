// -----------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2020, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2020, Knut Reinert & MPI für molekulare Genetik
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at: https://github.com/seqan/seqan3/blob/master/LICENSE.md
// -----------------------------------------------------------------------------------------------------

/*!\file
 * \brief Meta-header for the parallel module.
 * \author Rene Rahn <rene.rahn AT fu-berlin.de>
 */

#pragma once

/*!\defgroup parallel Parallel
 * \brief This module contains types and utilities for concurrent execution of algorithms in SeqAn.
 * \ingroup core
 *
 * \details
 *
 * ### Execution policies
 *
 * Since C++17/20 the standard defines execution policies for sequential, parallel and vectorised execution of
 * algorithms. These policies are not yet fully supported by the gcc compilers or only with a hard dependency on
 * external libraries. Thus, we define our own execution policies which adopt the behaviour of the policies defined
 * in the standard library. Once they are fully supported by the compilers our policies will merely alias the standard
 * policies.
 *
 * \if DEV
 * ### Concurrency support
 *
 * This module contains helper classes to synchronise threads in concurrent environments.
 *
 * \endif
 */

#include <seqan3/core/parallel/execution.hpp>
#include <seqan3/core/parallel/detail/all.hpp>
