// -----------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2023, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2023, Knut Reinert & MPI für molekulare Genetik
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at: https://github.com/seqan/seqan3/blob/master/LICENSE.md
// -----------------------------------------------------------------------------------------------------

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
