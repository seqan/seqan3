// SPDX-FileCopyrightText: 2006-2025 Knut Reinert & Freie Universität Berlin
// SPDX-FileCopyrightText: 2016-2025 Knut Reinert & MPI für molekulare Genetik
// SPDX-License-Identifier: BSD-3-Clause

/*!\file
 * \brief Provides a test literal "AC-GT"_gapped_dna5 for std::vector<seqan3::gapped<seqan3::dna5>>.
 * \author Marcel Ehrhardt <marcel.ehrhardt AT fu-berlin.de>
 */

#pragma once

#include <seqan3/alphabet/gap/gapped.hpp>
#include <seqan3/alphabet/nucleotide/dna5.hpp>

namespace seqan3::test
{

SEQAN3_WORKAROUND_LITERAL std::vector<seqan3::gapped<seqan3::dna5>> operator""_gapped_dna5(char const * s,
                                                                                           std::size_t n)
{
    std::vector<seqan3::gapped<seqan3::dna5>> r;
    r.resize(n);

    for (size_t i = 0; i < n; ++i)
        r[i].assign_char(s[i]);

    return r;
}

} // namespace seqan3::test
