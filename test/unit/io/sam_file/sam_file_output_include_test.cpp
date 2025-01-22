// SPDX-FileCopyrightText: 2006-2025 Knut Reinert & Freie Universität Berlin
// SPDX-FileCopyrightText: 2016-2025 Knut Reinert & MPI für molekulare Genetik
// SPDX-License-Identifier: BSD-3-Clause

#include <gtest/gtest.h>

// seqan3/io/sam_file/format_bam.hpp did not include core/debug_stream/tuple.hpp
// having only the following include leads to a compile error
#include <seqan3/io/sam_file/output.hpp>

TEST(sam_file_output, include)
{
    using sam_file_output_t = seqan3::sam_file_output<seqan3::fields<seqan3::field::id>,
                                                      seqan3::type_list<seqan3::format_sam, seqan3::format_bam>,
                                                      std::vector<std::string>>;
    std::string buffer{};
    std::ostringstream stream{buffer};
    {
        sam_file_output_t out{stream, std::vector<std::string>{}, std::vector<size_t>{}, seqan3::format_sam{}};
        out.emplace_back(std::string{});
    }
    {
        sam_file_output_t out{stream, std::vector<std::string>{}, std::vector<size_t>{}, seqan3::format_bam{}};
        out.emplace_back(std::string{});
    }
}
