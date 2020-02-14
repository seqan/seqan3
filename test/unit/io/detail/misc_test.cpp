// -----------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2020, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2020, Knut Reinert & MPI für molekulare Genetik
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at: https://github.com/seqan/seqan3/blob/master/LICENSE.md
// -----------------------------------------------------------------------------------------------------

#include <gtest/gtest.h>

#include <string>
#include <vector>

#include <seqan3/io/detail/magic_header.hpp>
#include <seqan3/io/detail/misc.hpp>
#include <seqan3/std/ranges>
#include <seqan3/test/tmp_filename.hpp>

struct dummy_file
{

    struct format1
    {
        static inline std::vector<std::string> file_extensions{ {"fa"}, {"fasta"}};
    };

    struct format2
    {
        static inline std::vector<std::string> file_extensions{ {"sam"}, {"bam"}};
    };

    using valid_formats = seqan3::type_list<format1, format2>;
};

TEST(misc, valid_file_extensions)
{
    // get all extensions.
    auto all_extensions = seqan3::detail::valid_file_extensions<dummy_file::valid_formats>();

    // define testing lambda
    auto cmp_lambda = [&all_extensions] (auto & source)
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

TEST(misc, valid_compression_extensions)
{
    using string_vector = std::vector<std::string>;
    string_vector valid_compression = seqan3::detail::valid_file_extensions<seqan3::detail::compression_formats>();

    #if defined(SEQAN3_HAS_ZLIB)
    // expect gz and bgzf
        EXPECT_TRUE(std::find(valid_compression.begin(), valid_compression.end(), "gz") != valid_compression.end());
        EXPECT_TRUE(std::find(valid_compression.begin(), valid_compression.end(), "bgzf") != valid_compression.end());
    #endif

    #if defined(SEQAN3_HAS_BZIP2)
        EXPECT_TRUE(std::find(valid_compression.begin(), valid_compression.end(), "bz2") != valid_compression.end());
    #endif

    #if defined(SEQAN3_HAS_ZSTD)
        EXPECT_TRUE(std::find(valid_compression.begin(), valid_compression.end(), "zst") != valid_compression.end());
    #endif
}
