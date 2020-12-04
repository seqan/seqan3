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
 * \ingroup utility
 *
 * \details
 *
 * ### Execution policies
 *
 * Here are currently only implementations which are part of detail and therefore not of interest for the common user. 
 *
 * \if DEV
 * ### Concurrency support
 *
 * This module contains helper classes to synchronise threads in concurrent environments.
 *
 * \endif
 */

 #include <seqan3/utility/parallel/detail/latch.hpp>
 #include <seqan3/utility/parallel/detail/reader_writer_manager.hpp>
 #include <seqan3/utility/parallel/detail/spin_delay.hpp>
