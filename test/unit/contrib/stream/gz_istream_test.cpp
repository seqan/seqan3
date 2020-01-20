// -----------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2020, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2020, Knut Reinert & MPI für molekulare Genetik
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at: https://github.com/seqan/seqan3/blob/master/LICENSE.md
// -----------------------------------------------------------------------------------------------------

#include <gtest/gtest.h>

#include <seqan3/contrib/stream/gz_istream.hpp>

#include "../../io/stream/istream_test_template.hpp"

using namespace seqan3;

template <>
class istream<contrib::gz_istream> : public ::testing::Test
{
public:
    static inline std::string compressed
    {
        '\x1f','\x8b','\x08','\x00','\x00','\x00','\x00','\x00','\x00','\x03','\x0b','\xc9','\x48','\x55','\x28','\x2c',
        '\xcd','\x4c','\xce','\x56','\x48','\x2a','\xca','\x2f','\xcf','\x53','\x48','\xcb','\xaf','\x50','\xc8','\x2a',
        '\xcd','\x2d','\x28','\x56','\xc8','\x2f','\x4b','\x2d','\x52','\x28','\x01','\x4a','\xe7','\x24','\x56','\x55',
        '\x2a','\xa4','\xe4','\xa7','\x03','\x00','\x39','\xa3','\x4f','\x41','\x2b','\x00','\x00','\x00'
    };
};

using test_types = ::testing::Types<contrib::gz_istream>;

INSTANTIATE_TYPED_TEST_SUITE_P(contrib_streams, istream, test_types, );
