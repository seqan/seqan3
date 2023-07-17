// -----------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2023, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2023, Knut Reinert & MPI für molekulare Genetik
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at: https://github.com/seqan/seqan3/blob/master/LICENSE.md
// -----------------------------------------------------------------------------------------------------

/*!\file
 * \brief Provides a test literal "1S1M1D1M1I"_cigar for std::vector<seqan3::cigar>.
 * \author Marcel Ehrhardt <marcel.ehrhardt AT fu-berlin.de>
 */

#pragma once

#include <string_view>

#include <seqan3/io/sam_file/detail/cigar.hpp>

namespace seqan3::test
{

SEQAN3_WORKAROUND_LITERAL std::vector<cigar> operator""_cigar(char const * s, std::size_t n)
{
    return seqan3::detail::parse_cigar(std::string_view{s, n});
}

} // namespace seqan3::test
