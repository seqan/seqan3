// -----------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2019, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2019, Knut Reinert & MPI für molekulare Genetik
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at: https://github.com/seqan/seqan3/blob/master/LICENSE.md
// -----------------------------------------------------------------------------------------------------

#include <gtest/gtest.h>

#include <meta/meta.hpp>

#include <seqan3/alphabet/aminoacid/all.hpp>
#include <seqan3/core/type_list.hpp>

using namespace seqan3;

template <typename T>
using aminoacid_conversion = ::testing::Test;

using aminoacid_types = type_list<aa10li, aa10murphy, aa20, aa27>; // needed for some tests
using aminoacid_gtest_types = detail::transfer_template_args_onto_t<aminoacid_types, ::testing::Types>;

TYPED_TEST_CASE(aminoacid_conversion, aminoacid_gtest_types);

// conversion to any other amino acid type
TYPED_TEST(aminoacid_conversion, explicit_conversion)
{
    meta::for_each(aminoacid_types{}, [&] (auto && aa) constexpr
    {
        using out_type = std::decay_t<decltype(aa)>;
        EXPECT_EQ(static_cast<out_type>(TypeParam{}.assign_char('A')), out_type{}.assign_char('A'));
        EXPECT_EQ(static_cast<out_type>(TypeParam{}.assign_char('C')), out_type{}.assign_char('C'));
        EXPECT_EQ(static_cast<out_type>(TypeParam{}.assign_char('F')), out_type{}.assign_char('F'));
        EXPECT_EQ(static_cast<out_type>(TypeParam{}.assign_char('G')), out_type{}.assign_char('G'));
        EXPECT_EQ(static_cast<out_type>(TypeParam{}.assign_char('H')), out_type{}.assign_char('H'));
        EXPECT_EQ(static_cast<out_type>(TypeParam{}.assign_char('I')), out_type{}.assign_char('I'));
        EXPECT_EQ(static_cast<out_type>(TypeParam{}.assign_char('K')), out_type{}.assign_char('K'));
        EXPECT_EQ(static_cast<out_type>(TypeParam{}.assign_char('P')), out_type{}.assign_char('P'));

        if (std::is_same_v<TypeParam, aa27>)
        {
            EXPECT_EQ(static_cast<out_type>(TypeParam{}.assign_char('D')), out_type{}.assign_char('D'));
            EXPECT_EQ(static_cast<out_type>(TypeParam{}.assign_char('E')), out_type{}.assign_char('E'));
            EXPECT_EQ(static_cast<out_type>(TypeParam{}.assign_char('L')), out_type{}.assign_char('L'));
            EXPECT_EQ(static_cast<out_type>(TypeParam{}.assign_char('M')), out_type{}.assign_char('M'));
            EXPECT_EQ(static_cast<out_type>(TypeParam{}.assign_char('N')), out_type{}.assign_char('N'));
            EXPECT_EQ(static_cast<out_type>(TypeParam{}.assign_char('Q')), out_type{}.assign_char('Q'));
            EXPECT_EQ(static_cast<out_type>(TypeParam{}.assign_char('R')), out_type{}.assign_char('R'));
            EXPECT_EQ(static_cast<out_type>(TypeParam{}.assign_char('S')), out_type{}.assign_char('S'));
            EXPECT_EQ(static_cast<out_type>(TypeParam{}.assign_char('T')), out_type{}.assign_char('T'));
            EXPECT_EQ(static_cast<out_type>(TypeParam{}.assign_char('V')), out_type{}.assign_char('V'));
            EXPECT_EQ(static_cast<out_type>(TypeParam{}.assign_char('W')), out_type{}.assign_char('W'));
            EXPECT_EQ(static_cast<out_type>(TypeParam{}.assign_char('Y')), out_type{}.assign_char('Y'));
            EXPECT_EQ(static_cast<out_type>(TypeParam{}.assign_char('B')), out_type{}.assign_char('B'));
            EXPECT_EQ(static_cast<out_type>(TypeParam{}.assign_char('J')), out_type{}.assign_char('J'));
            EXPECT_EQ(static_cast<out_type>(TypeParam{}.assign_char('O')), out_type{}.assign_char('O'));
            EXPECT_EQ(static_cast<out_type>(TypeParam{}.assign_char('U')), out_type{}.assign_char('U'));
            EXPECT_EQ(static_cast<out_type>(TypeParam{}.assign_char('X')), out_type{}.assign_char('X'));
            EXPECT_EQ(static_cast<out_type>(TypeParam{}.assign_char('Z')), out_type{}.assign_char('Z'));
            EXPECT_EQ(static_cast<out_type>(TypeParam{}.assign_char('*')), out_type{}.assign_char('*'));
            EXPECT_EQ(static_cast<out_type>(TypeParam{}.assign_char('!')), out_type{}.assign_char('!'));
        }
        else if (std::is_same_v<TypeParam, aa20>)
        {
            EXPECT_EQ(static_cast<out_type>(TypeParam{}.assign_char('D')), out_type{}.assign_char('D'));
            EXPECT_EQ(static_cast<out_type>(TypeParam{}.assign_char('E')), out_type{}.assign_char('E'));
            EXPECT_EQ(static_cast<out_type>(TypeParam{}.assign_char('L')), out_type{}.assign_char('L'));
            EXPECT_EQ(static_cast<out_type>(TypeParam{}.assign_char('M')), out_type{}.assign_char('M'));
            EXPECT_EQ(static_cast<out_type>(TypeParam{}.assign_char('N')), out_type{}.assign_char('N'));
            EXPECT_EQ(static_cast<out_type>(TypeParam{}.assign_char('Q')), out_type{}.assign_char('Q'));
            EXPECT_EQ(static_cast<out_type>(TypeParam{}.assign_char('R')), out_type{}.assign_char('R'));
            EXPECT_EQ(static_cast<out_type>(TypeParam{}.assign_char('S')), out_type{}.assign_char('S'));
            EXPECT_EQ(static_cast<out_type>(TypeParam{}.assign_char('T')), out_type{}.assign_char('T'));
            EXPECT_EQ(static_cast<out_type>(TypeParam{}.assign_char('V')), out_type{}.assign_char('V'));
            EXPECT_EQ(static_cast<out_type>(TypeParam{}.assign_char('W')), out_type{}.assign_char('W'));
            EXPECT_EQ(static_cast<out_type>(TypeParam{}.assign_char('Y')), out_type{}.assign_char('Y'));
            EXPECT_EQ(static_cast<out_type>(TypeParam{}.assign_char('B')), out_type{}.assign_char('D'));
            EXPECT_EQ(static_cast<out_type>(TypeParam{}.assign_char('J')), out_type{}.assign_char('L'));
            EXPECT_EQ(static_cast<out_type>(TypeParam{}.assign_char('O')), out_type{}.assign_char('L'));
            EXPECT_EQ(static_cast<out_type>(TypeParam{}.assign_char('U')), out_type{}.assign_char('C'));
            EXPECT_EQ(static_cast<out_type>(TypeParam{}.assign_char('X')), out_type{}.assign_char('S'));
            EXPECT_EQ(static_cast<out_type>(TypeParam{}.assign_char('Z')), out_type{}.assign_char('E'));
            EXPECT_EQ(static_cast<out_type>(TypeParam{}.assign_char('*')), out_type{}.assign_char('W'));
            EXPECT_EQ(static_cast<out_type>(TypeParam{}.assign_char('!')), out_type{}.assign_char('S'));
        }
        else if (std::is_same_v<TypeParam, aa10murphy>)
        {
            EXPECT_EQ(static_cast<out_type>(TypeParam{}.assign_char('D')), out_type{}.assign_char('B'));
            EXPECT_EQ(static_cast<out_type>(TypeParam{}.assign_char('E')), out_type{}.assign_char('B'));
            EXPECT_EQ(static_cast<out_type>(TypeParam{}.assign_char('J')), out_type{}.assign_char('I'));
            EXPECT_EQ(static_cast<out_type>(TypeParam{}.assign_char('L')), out_type{}.assign_char('I'));
            EXPECT_EQ(static_cast<out_type>(TypeParam{}.assign_char('M')), out_type{}.assign_char('I'));
            EXPECT_EQ(static_cast<out_type>(TypeParam{}.assign_char('N')), out_type{}.assign_char('B'));
            EXPECT_EQ(static_cast<out_type>(TypeParam{}.assign_char('O')), out_type{}.assign_char('K'));
            EXPECT_EQ(static_cast<out_type>(TypeParam{}.assign_char('Q')), out_type{}.assign_char('B'));
            EXPECT_EQ(static_cast<out_type>(TypeParam{}.assign_char('R')), out_type{}.assign_char('K'));
            EXPECT_EQ(static_cast<out_type>(TypeParam{}.assign_char('S')), out_type{}.assign_char('S'));
            EXPECT_EQ(static_cast<out_type>(TypeParam{}.assign_char('T')), out_type{}.assign_char('S'));
            EXPECT_EQ(static_cast<out_type>(TypeParam{}.assign_char('U')), out_type{}.assign_char('C'));
            EXPECT_EQ(static_cast<out_type>(TypeParam{}.assign_char('V')), out_type{}.assign_char('I'));
            EXPECT_EQ(static_cast<out_type>(TypeParam{}.assign_char('W')), out_type{}.assign_char('F'));
            EXPECT_EQ(static_cast<out_type>(TypeParam{}.assign_char('X')), out_type{}.assign_char('S'));
            EXPECT_EQ(static_cast<out_type>(TypeParam{}.assign_char('Y')), out_type{}.assign_char('F'));
            EXPECT_EQ(static_cast<out_type>(TypeParam{}.assign_char('Z')), out_type{}.assign_char('B'));
            EXPECT_EQ(static_cast<out_type>(TypeParam{}.assign_char('*')), out_type{}.assign_char('F'));
            EXPECT_EQ(static_cast<out_type>(TypeParam{}.assign_char('!')), out_type{}.assign_char('S'));
        }
        else if (std::is_same_v<TypeParam, aa10li>)
        {
            EXPECT_EQ(static_cast<out_type>(TypeParam{}.assign_char('D')), out_type{}.assign_char('B'));
            EXPECT_EQ(static_cast<out_type>(TypeParam{}.assign_char('E')), out_type{}.assign_char('B'));
            EXPECT_EQ(static_cast<out_type>(TypeParam{}.assign_char('J')), out_type{}.assign_char('J'));
            EXPECT_EQ(static_cast<out_type>(TypeParam{}.assign_char('L')), out_type{}.assign_char('J'));
            EXPECT_EQ(static_cast<out_type>(TypeParam{}.assign_char('M')), out_type{}.assign_char('J'));
            EXPECT_EQ(static_cast<out_type>(TypeParam{}.assign_char('N')), out_type{}.assign_char('H'));
            EXPECT_EQ(static_cast<out_type>(TypeParam{}.assign_char('O')), out_type{}.assign_char('K'));
            EXPECT_EQ(static_cast<out_type>(TypeParam{}.assign_char('Q')), out_type{}.assign_char('B'));
            EXPECT_EQ(static_cast<out_type>(TypeParam{}.assign_char('R')), out_type{}.assign_char('K'));
            EXPECT_EQ(static_cast<out_type>(TypeParam{}.assign_char('S')), out_type{}.assign_char('A'));
            EXPECT_EQ(static_cast<out_type>(TypeParam{}.assign_char('T')), out_type{}.assign_char('A'));
            EXPECT_EQ(static_cast<out_type>(TypeParam{}.assign_char('U')), out_type{}.assign_char('C'));
            EXPECT_EQ(static_cast<out_type>(TypeParam{}.assign_char('V')), out_type{}.assign_char('I'));
            EXPECT_EQ(static_cast<out_type>(TypeParam{}.assign_char('W')), out_type{}.assign_char('F'));
            EXPECT_EQ(static_cast<out_type>(TypeParam{}.assign_char('X')), out_type{}.assign_char('A'));
            EXPECT_EQ(static_cast<out_type>(TypeParam{}.assign_char('Y')), out_type{}.assign_char('F'));
            EXPECT_EQ(static_cast<out_type>(TypeParam{}.assign_char('Z')), out_type{}.assign_char('B'));
            EXPECT_EQ(static_cast<out_type>(TypeParam{}.assign_char('*')), out_type{}.assign_char('F'));
            EXPECT_EQ(static_cast<out_type>(TypeParam{}.assign_char('!')), out_type{}.assign_char('A'));
        }
    });
}
