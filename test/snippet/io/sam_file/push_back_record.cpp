// SPDX-FileCopyrightText: 2006-2025 Knut Reinert & Freie Universität Berlin
// SPDX-FileCopyrightText: 2016-2025 Knut Reinert & MPI für molekulare Genetik
// SPDX-License-Identifier: CC0-1.0

#include <sstream>
#include <string>

#include <seqan3/io/sam_file/output.hpp>
#include <seqan3/utility/type_list/type_list.hpp>

int main()
{
    seqan3::sam_file_output fout{std::ostringstream{}, seqan3::format_sam{}};

    seqan3::record<seqan3::type_list<uint32_t, std::string>, seqan3::fields<seqan3::field::mapq, seqan3::field::id>> r;

    // ...

    fout.push_back(r);
}
