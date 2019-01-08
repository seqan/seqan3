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
 * \brief Provides pairwise alignment function.
 * \author Rene Rahn <rene.rahn AT fu-berlin.de>
 */

#pragma once

#include <functional>
#include <iostream>
#include <tuple>
#include <type_traits>

#include <meta/meta.hpp>

#include <range/v3/view/bounded.hpp>
#include <range/v3/view/transform.hpp>
#include <range/v3/view/single.hpp>

#include <seqan3/alignment/configuration/all.hpp>
#include <seqan3/alignment/pairwise/align_result.hpp>
#include <seqan3/alignment/pairwise/alignment_configurator.hpp>
#include <seqan3/alignment/pairwise/execution/all.hpp>

#include <seqan3/alphabet/gap/gapped.hpp>

#include <seqan3/core/algorithm/all.hpp>
#include <seqan3/core/metafunction/basic.hpp>
#include <seqan3/range/view/persist.hpp>

#include <seqan3/std/concepts>
#include <seqan3/std/ranges>

namespace seqan3
{

//!\cond
template <std::ranges::InputRange sequence_t, typename alignment_config_t>
    requires detail::is_type_specialisation_of_v<remove_cvref_t<alignment_config_t>, configuration> &&
             tuple_like_concept<value_type_t<std::ranges::iterator_t<std::remove_reference_t<sequence_t>>>>
constexpr auto align_pairwise(sequence_t && seq, alignment_config_t const & config)
{
    // Wrap in persist to make the code also work with temporaries that are not views.
    auto seq_view = view::persist(std::forward<sequence_t>(seq));
    // Configure the alignment algorithm.
    auto kernel = detail::alignment_configurator::configure(seq_view, config);
    // Create a two-way executor for the alignment.
    detail::alignment_executor_two_way exec{std::move(seq_view), kernel};
    // Return the range over the alignments.
    return alignment_range{std::move(exec)};
}

template <tuple_like_concept seq_t,
          typename alignment_config_t>
    requires detail::is_type_specialisation_of_v<remove_cvref_t<alignment_config_t>, configuration>
constexpr auto align_pairwise(seq_t && seq, alignment_config_t && config)
{
    static_assert(std::tuple_size_v<std::remove_reference_t<seq_t>> == 2,
                  "Expects exactly two sequences for pairwise alignments.");

    return align_pairwise(ranges::view::single(std::forward<seq_t>(seq)) | ranges::view::bounded,
                          std::forward<alignment_config_t>(config));
}
//!\endcond
} // namespace seqan3
