// SPDX-FileCopyrightText: 2006-2025 Knut Reinert & Freie Universität Berlin
// SPDX-FileCopyrightText: 2016-2025 Knut Reinert & MPI für molekulare Genetik
// SPDX-License-Identifier: BSD-3-Clause

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
