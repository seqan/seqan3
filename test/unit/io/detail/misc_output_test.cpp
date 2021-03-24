// -----------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2020, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2020, Knut Reinert & MPI für molekulare Genetik
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at: https://github.com/seqan/seqan3/blob/master/LICENSE.md
// -----------------------------------------------------------------------------------------------------

#include <gtest/gtest.h>

#include <fstream>
#include <vector>

#include <seqan3/core/debug_stream.hpp>

#include <seqan3/io/detail/misc_input.hpp>
#include <seqan3/io/detail/misc_output.hpp>
#include <seqan3/test/tmp_filename.hpp>

inline seqan3::test::tmp_filename tmp_compressed_file(std::string const & file_extension)
{
    std::string const tmp_file_name = "io_misc_output_test.txt." + file_extension;
    seqan3::test::tmp_filename tmp_file{tmp_file_name.c_str()};

    // We need a copy of the path because `make_secondary_ostream` will strip the compression extension.
    auto file_path = tmp_file.get_path();
    std::ofstream filestream{file_path};
    auto stream_ptr = seqan3::detail::make_secondary_ostream(filestream, file_path);
    *stream_ptr << std::string(8, 'a') << '\n';

    return tmp_file;
}

inline std::vector<char> read_file_content(std::filesystem::path const & path)
{
    std::ifstream filestream{path};
    using char_t = decltype(filestream)::char_type;
    return {std::istreambuf_iterator{filestream}, std::istreambuf_iterator<char_t>{}};
}

#ifdef SEQAN3_HAS_ZLIB
TEST(misc_output, issue2455_gz)
{
    seqan3::test::tmp_filename const compressed_file = tmp_compressed_file("gz");
    std::vector<char> const file_content = read_file_content(compressed_file.get_path());

    EXPECT_TRUE(seqan3::detail::starts_with(file_content, seqan3::detail::gz_compression::magic_header));
    // gz should not have a valid bgzf header (the gz header is a prefix of the bgzf header)
    EXPECT_FALSE(seqan3::detail::bgzf_compression::validate_header(std::span{file_content}));
}

TEST(misc_output, issue2455_bgzf)
{
    seqan3::test::tmp_filename const compressed_file = tmp_compressed_file("bgzf");
    std::vector<char> const file_content = read_file_content(compressed_file.get_path());

    EXPECT_TRUE(seqan3::detail::bgzf_compression::validate_header(std::span{file_content}));
}
#endif

#ifdef SEQAN3_HAS_BZIP2
TEST(misc_output, issue2455_bz)
{
    seqan3::test::tmp_filename const compressed_file = tmp_compressed_file("bz2");
    std::vector<char> const file_content = read_file_content(compressed_file.get_path());

    EXPECT_TRUE(seqan3::detail::starts_with(file_content, seqan3::detail::bz2_compression::magic_header));
}
#endif
