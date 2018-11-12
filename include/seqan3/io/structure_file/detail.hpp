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
 * \brief Helper functions (e.g. conversions) for the structure IO submodule.
 * \author Jörg Winkler <j.winkler AT fu-berlin.de>
 */

#pragma once

#include <map>
#include <stack>

#include <seqan3/alphabet/structure/RnaStructure.hpp>

namespace seqan3::detail
{
/*!
 * \brief Transforms a structure annotation string into a base pair probability matrix.
 * \ingroup structure_file
 * \throws seqan3::parse_error if unpaired brackets are found in the structure annotation.
 * \tparam structure_alph_type The type of the structure alphabet; must satisfy seqan3::RnaStructure.
 * \tparam bpp_type            The type of the target matrix.
 * \tparam structure_type      The range type of the structure annotation.
 * \param[out] bpp             The target matrix that receives the base pair probabilities.
 * \param[in]  structure       The source structure annotation.
 * \param[in]  weight          The weight to be assigned to all interactions present.
 *                             As the source allows only one interaction partner, the weight defaults to 1.0.
 */
template <typename structure_alph_type, typename bpp_type, std::ranges::Range structure_type>
inline
void bpp_from_rna_structure(bpp_type & bpp, structure_type const & structure, double weight = 1.)
{
    if constexpr (!RnaStructure<structure_alph_type>)
        throw parse_error{"Cannot create base pair probabilities from a structure that is not RNA structure."};

    bpp.clear();
    if constexpr (std::ranges::SizedRange<structure_type>)
        bpp.reserve(size(structure));

    std::stack<size_t> brackets[max_pseudoknot_depth_v<structure_alph_type>];
    size_t pos = 0ul;
    for (structure_alph_type symbol : structure)
    {
        bpp.push_back({});
        uint8_t const id = pseudoknot_id(symbol).value_or(0);

        if (symbol.is_pair_open())
        {
            brackets[id].push(pos);
        }
        else if (symbol.is_pair_close())
        {
            if (!brackets[id].empty())
            {
                bpp[pos].emplace(weight, brackets[id].top());
                bpp[brackets[id].top()].emplace(weight, pos);
                brackets[id].pop();
            }
            else
            {
                throw parse_error{std::string{"Invalid bracket notation: Unpaired closing bracket at position "}
                                  + std::to_string(pos) + "."};
            };
        }
        // no actions for unpaired
        ++pos;
    }
    for (uint8_t id = 0u; id < max_pseudoknot_depth_v<structure_alph_type>; ++id)
    {
        if (!brackets[id].empty())
        {
            throw parse_error{std::string{"Invalid bracket notation: Unpaired opening bracket at position "}
                              + std::to_string(brackets[id].top()) + "."};
        }
    }
}

} // namespace seqan3::detail
