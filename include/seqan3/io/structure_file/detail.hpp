// -----------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2020, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2020, Knut Reinert & MPI für molekulare Genetik
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at: https://github.com/seqan/seqan3/blob/master/LICENSE.md
// -----------------------------------------------------------------------------------------------------

/*!\file
 * \brief Helper functions (e.g. conversions) for the structure IO submodule.
 * \author Jörg Winkler <j.winkler AT fu-berlin.de>
 */

#pragma once

#include <map>
#include <stack>

#include <seqan3/alphabet/structure/concept.hpp>
#include <seqan3/io/exception.hpp>
#include <seqan3/std/ranges>

namespace seqan3::detail
{
/*!
 * \brief Transforms a structure annotation string into a base pair probability matrix.
 * \ingroup structure_file
 * \throws seqan3::parse_error if unpaired brackets are found in the structure annotation.
 * \tparam structure_alph_type The type of the structure alphabet; must satisfy seqan3::rna_structure_alphabet.
 * \tparam bpp_type            The type of the target matrix.
 * \tparam structure_type      The range type of the structure annotation.
 * \param[out] bpp             The target matrix that receives the base pair probabilities.
 * \param[in]  structure       The source structure annotation.
 * \param[in]  weight          The weight to be assigned to all interactions present.
 *                             As the source allows only one interaction partner, the weight defaults to 1.0.
 */
template <typename structure_alph_type, typename bpp_type, std::ranges::range structure_type>
inline
void bpp_from_rna_structure(bpp_type & bpp, structure_type const & structure, double weight = 1.)
{
    if constexpr (!rna_structure_alphabet<structure_alph_type>)
        throw parse_error{"Cannot create base pair probabilities from a structure that is not RNA structure."};

    bpp.clear();
    if constexpr (std::ranges::sized_range<structure_type>)
        bpp.reserve(size(structure));

    std::stack<size_t> brackets[max_pseudoknot_depth<structure_alph_type>];
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
    for (uint8_t id = 0u; id < max_pseudoknot_depth<structure_alph_type>; ++id)
    {
        if (!brackets[id].empty())
        {
            throw parse_error{std::string{"Invalid bracket notation: Unpaired opening bracket at position "}
                              + std::to_string(brackets[id].top()) + "."};
        }
    }
}

} // namespace seqan3::detail
