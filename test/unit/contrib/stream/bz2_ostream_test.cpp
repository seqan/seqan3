// -----------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2020, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2020, Knut Reinert & MPI für molekulare Genetik
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at: https://github.com/seqan/seqan3/blob/master/LICENSE.md
// -----------------------------------------------------------------------------------------------------

#include <gtest/gtest.h>

#include <seqan3/contrib/stream/bz2_ostream.hpp>

#include "../../io/stream/ostream_test_template.hpp"

template <>
class ostream<seqan3::contrib::bz2_ostream> : public ::testing::Test
{
public:
    static inline std::string compressed
    {
        '\x42','\x5A','\x68','\x39','\x31','\x41','\x59','\x26','\x53','\x59','\x45','\x9D','\xEE','\x61','\x00','\x00',
        '\x04','\x13','\x80','\x40','\x00','\x04','\x00','\x3F','\xFF','\xFF','\xF0','\x20','\x00','\x31','\x46','\x86',
        '\x80','\x00','\x00','\x31','\xE9','\xA9','\xA6','\x4C','\x86','\x11','\xB4','\x6D','\x47','\x62','\x62','\x08',
        '\x49','\xED','\x7A','\xA1','\x53','\x65','\x65','\xB1','\x25','\xE3','\xE2','\x60','\xB1','\xF8','\x98','\x39',
        '\xDD','\x4C','\x09','\x6F','\x9C','\xE8','\x5D','\xC9','\x14','\xE1','\x42','\x41','\x16','\x77','\xB9','\x84',
    };
};

using test_types = ::testing::Types<seqan3::contrib::bz2_ostream>;

INSTANTIATE_TYPED_TEST_SUITE_P(contrib_streams, ostream, test_types, );
