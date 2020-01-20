// -----------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2020, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2020, Knut Reinert & MPI für molekulare Genetik
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at: https://github.com/seqan/seqan3/blob/master/LICENSE.md
// -----------------------------------------------------------------------------------------------------

#include <gtest/gtest.h>

#include <seqan3/contrib/stream/bgzf_istream.hpp>

#include "../../io/stream/istream_test_template.hpp"

using namespace seqan3;

template <>
class istream<contrib::bgzf_istream> : public ::testing::Test
{
public:
    static inline std::string compressed
    {
        '\x1F', '\x8B', '\x08', '\x04', '\x00', '\x00', '\x00', '\x00', '\x00', '\xFF', '\x06', '\x00', '\x42', '\x43',
        '\x02', '\x00', '\x45', '\x00', '\x0B', '\xC9', '\x48', '\x55', '\x28', '\x2C', '\xCD', '\x4C', '\xCE', '\x56',
        '\x48', '\x2A', '\xCA', '\x2F', '\xCF', '\x53', '\x48', '\xCB', '\xAF', '\x50', '\xC8', '\x2A', '\xCD', '\x2D',
        '\x28', '\x56', '\xC8', '\x2F', '\x4B', '\x2D', '\x52', '\x28', '\x01', '\x4A', '\xE7', '\x24', '\x56', '\x55',
        '\x2A', '\xA4', '\xE4', '\xA7', '\x03', '\x00', '\x39', '\xA3', '\x4F', '\x41', '\x2B', '\x00', '\x00', '\x00',
        '\x1F', '\x8B', '\x08', '\x04', '\x00', '\x00', '\x00', '\x00', '\x00', '\xFF', '\x06', '\x00', '\x42', '\x43',
        '\x02', '\x00', '\x1B', '\x00', '\x03', '\x00', '\x00', '\x00', '\x00', '\x00', '\x00', '\x00', '\x00', '\x00'
    };
};

using test_types = ::testing::Types<contrib::bgzf_istream>;

INSTANTIATE_TYPED_TEST_SUITE_P(contrib_streams, istream, test_types, );
