// -----------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2023, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2023, Knut Reinert & MPI für molekulare Genetik
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at: https://github.com/seqan/seqan3/blob/master/LICENSE.md
// -----------------------------------------------------------------------------------------------------

#include <gtest/gtest.h>

#include <string>
#include <vector>

#include <seqan3/io/detail/misc.hpp>

struct dummy_file
{

    struct format1
    {
        static inline std::vector<std::string> file_extensions{{"fa"}, {"fasta"}};
    };

    struct format2
    {
        static inline std::vector<std::string> file_extensions{{"sam"}, {"bam"}};
    };

    using valid_formats = seqan3::type_list<format1, format2>;
};

TEST(misc, valid_file_extensions)
{
    // get all extensions.
    auto all_extensions = seqan3::detail::valid_file_extensions<dummy_file::valid_formats>();

    // define testing lambda
    auto cmp_lambda = [&all_extensions](auto & source)
    {
        return std::find(all_extensions.begin(), all_extensions.end(), source);
    };

    // Test format1 extensions
    for (std::string & ext : dummy_file::format1::file_extensions)
        EXPECT_NE(cmp_lambda(ext), all_extensions.end());

    // Test format2 extensions
    for (std::string & ext : dummy_file::format2::file_extensions)
        EXPECT_NE(cmp_lambda(ext), all_extensions.end());
}
