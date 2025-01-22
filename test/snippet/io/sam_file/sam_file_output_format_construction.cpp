// SPDX-FileCopyrightText: 2006-2025 Knut Reinert & Freie Universität Berlin
// SPDX-FileCopyrightText: 2016-2025 Knut Reinert & MPI für molekulare Genetik
// SPDX-License-Identifier: CC0-1.0

#include <sstream>

#include <seqan3/io/sam_file/output.hpp>

int main()
{
    // no need to specify the template arguments <...> for format specialization:
    seqan3::sam_file_output fout{std::ostringstream{}, seqan3::format_sam{}, seqan3::fields<seqan3::field::mapq>{}};
}
