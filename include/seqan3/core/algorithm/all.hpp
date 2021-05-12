// -----------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2020, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2020, Knut Reinert & MPI für molekulare Genetik
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at: https://github.com/seqan/seqan3/blob/master/LICENSE.md
// -----------------------------------------------------------------------------------------------------

/*!\file
 * \brief Meta-Header for components of the algorithm submodule.
 * \author Rene Rahn <rene.rahn AT fu-berlin.de>
 */

#pragma once

/*!\defgroup algorithm Algorithm
 * \ingroup core
 * \brief Provides core functionality used to configure algorithms.
 */
#include <seqan3/core/algorithm/algorithm_result_generator_range.hpp>
#include <seqan3/core/algorithm/bound.hpp>
#if SEQAN3_VERSION_MAJOR == 3 && SEQAN3_VERSION_MINOR == 1
  #pragma warning "Remove #include <seqan3/core/algorithm/bound.hpp> from this header."
#endif
#include <seqan3/core/configuration/all.hpp>
#if SEQAN3_VERSION_MAJOR == 3 && SEQAN3_VERSION_MINOR == 1
  #pragma warning "Remove #include <seqan3/core/configuration/all.hpp> from this header."
#endif
