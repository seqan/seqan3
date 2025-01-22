// SPDX-FileCopyrightText: 2006-2025 Knut Reinert & Freie Universität Berlin
// SPDX-FileCopyrightText: 2016-2025 Knut Reinert & MPI für molekulare Genetik
// SPDX-License-Identifier: BSD-3-Clause

#include <gtest/gtest.h>

#include <seqan3/contrib/stream/bgzf_ostream.hpp>

#include "../../io/stream/ostream_test_template.hpp"

template <>
class ostream<seqan3::contrib::bgzf_ostream> : public ::testing::Test
{
public:
    static constexpr bool zero_out_os_byte = true;

    static inline std::string compressed{
        // OS = 0
        '\x1F', '\x8B', '\x08', '\x04', '\x00', '\x00', '\x00', '\x00', '\x00', '\x00', '\x06', '\x00', '\x42',
        '\x43', '\x02', '\x00', '\x45', '\x00', '\x0B', '\xC9', '\x48', '\x55', '\x28', '\x2C', '\xCD', '\x4C',
        '\xCE', '\x56', '\x48', '\x2A', '\xCA', '\x2F', '\xCF', '\x53', '\x48', '\xCB', '\xAF', '\x50', '\xC8',
        '\x2A', '\xCD', '\x2D', '\x28', '\x56', '\xC8', '\x2F', '\x4B', '\x2D', '\x52', '\x28', '\x01', '\x4A',
        '\xE7', '\x24', '\x56', '\x55', '\x2A', '\xA4', '\xE4', '\xA7', '\x03', '\x00', '\x39', '\xA3', '\x4F',
        '\x41', '\x2B', '\x00', '\x00', '\x00', '\x1F', '\x8B', '\x08', '\x04', '\x00', '\x00', '\x00', '\x00',
        '\x00', '\xFF', '\x06', '\x00', '\x42', '\x43', '\x02', '\x00', '\x1B', '\x00', '\x03', '\x00', '\x00',
        '\x00', '\x00', '\x00', '\x00', '\x00', '\x00', '\x00'}; // Note we zeroed the 10th byte which indicates the OS on which the file was compressed.
};

using test_types = ::testing::Types<seqan3::contrib::bgzf_ostream>;

INSTANTIATE_TYPED_TEST_SUITE_P(contrib_streams, ostream, test_types, );
