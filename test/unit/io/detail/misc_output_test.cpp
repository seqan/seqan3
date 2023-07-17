// -----------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2023, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2023, Knut Reinert & MPI für molekulare Genetik
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at: https://github.com/seqan/seqan3/blob/master/LICENSE.md
// -----------------------------------------------------------------------------------------------------

#include <gtest/gtest.h>

#include <fstream>
#include <vector>

#include <seqan3/core/debug_stream.hpp>
#include <seqan3/io/detail/misc_input.hpp>
#include <seqan3/io/detail/misc_output.hpp>
#include <seqan3/test/tmp_directory.hpp>

// We need a copy of the path because `make_secondary_ostream` will strip the compression extension.
inline void tmp_compressed_file(std::filesystem::path filename)
{
    std::ofstream filestream{filename};
    auto stream_ptr = seqan3::detail::make_secondary_ostream(filestream, filename);
    *stream_ptr << std::string(8, 'a') << '\n';
}

inline std::vector<char> read_file_content(std::filesystem::path const & path)
{
    std::ifstream filestream{path};
    using char_t = decltype(filestream)::char_type;
    return {std::istreambuf_iterator{filestream}, std::istreambuf_iterator<char_t>{}};
}

#if defined(SEQAN3_HAS_ZLIB)
TEST(misc_output, issue2455_gz)
{
    seqan3::test::tmp_directory tmp;
    auto filename = tmp.path() / "io_misc_output_test.txt.gz";
    tmp_compressed_file(filename);
    std::vector<char> const file_content = read_file_content(filename);

    EXPECT_TRUE(seqan3::detail::starts_with(file_content, seqan3::detail::gz_compression::magic_header));
    // gz should not have a valid bgzf header (the gz header is a prefix of the bgzf header)
    EXPECT_FALSE(seqan3::detail::bgzf_compression::validate_header(std::span{file_content}));
}

TEST(misc_output, issue2455_bgzf)
{
    seqan3::test::tmp_directory tmp;
    auto filename = tmp.path() / "io_misc_output_test.txt.bgzf";
    tmp_compressed_file(filename);
    std::vector<char> const file_content = read_file_content(filename);

    EXPECT_TRUE(seqan3::detail::bgzf_compression::validate_header(std::span{file_content}));
}
#endif

#if defined(SEQAN3_HAS_BZIP2)
TEST(misc_output, issue2455_bz)
{
    seqan3::test::tmp_directory tmp;
    auto filename = tmp.path() / "io_misc_output_test.txt.bz2";
    tmp_compressed_file(filename);
    std::vector<char> const file_content = read_file_content(filename);

    EXPECT_TRUE(seqan3::detail::starts_with(file_content, seqan3::detail::bz2_compression::magic_header));
}
#endif
