// ============================================================================
//                 SeqAn - The Library for Sequence Analysis
// ============================================================================
//
// Copyright (c) 2006-2018, Knut Reinert & Freie Universitaet Berlin
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
// ============================================================================

 /*!\file
  * \brief Meta-header for the \link configuration alignment configuration module \endlink.
  * \author Rene Rahn <rene.rahn AT fu-berlin.de>
  */

 #pragma once

/*!\defgroup configuration Configuration
 * \brief Data structures and utility functions for configuring alignment algorithm.
 * \ingroup alignment
 *
 * \todo Write detailed landing page.
 */

#include <seqan3/alignment/configuration/align_config_edit.hpp>
#include <seqan3/alignment/configuration/align_config_gap.hpp>
#include <seqan3/alignment/configuration/align_config_global.hpp>
#include <seqan3/alignment/configuration/align_config_max_error.hpp>
#include <seqan3/alignment/configuration/align_config_output.hpp>
#include <seqan3/alignment/configuration/align_config_score.hpp>
#include <seqan3/alignment/configuration/align_config_sequence_ends.hpp>
#include <seqan3/alignment/configuration/utility.hpp>

/*!\namespace seqan3::align_cfg
 * \brief A special sub namespace for the alignment configurations.
 */
