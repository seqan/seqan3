// SPDX-FileCopyrightText: 2006-2025 Knut Reinert & Freie Universität Berlin
// SPDX-FileCopyrightText: 2016-2025 Knut Reinert & MPI für molekulare Genetik
// SPDX-License-Identifier: CC0-1.0

#include <seqan3/core/debug_stream.hpp>
#include <seqan3/io/sam_file/sam_tag_dictionary.hpp>

int main()
{
    using namespace seqan3::literals;

    seqan3::sam_tag_dictionary dict{}; // initialise empty dictionary

    dict.get<"NM"_tag>() = 3;         // set SAM tag 'NM' to 3 (integer type)
    dict.get<"CO"_tag>() = "comment"; // set SAM tag 'CO' to "comment" (string type)

    auto nm = dict.get<"NM"_tag>(); // get SAM tag 'NM' (note: type is int32_t)
    auto co = dict.get<"CO"_tag>(); // get SAM tag 'CO' (note: type is std::string)

    seqan3::debug_stream << nm << '\n'; // will print '3'
    seqan3::debug_stream << co << '\n'; // will print "comment"
}
