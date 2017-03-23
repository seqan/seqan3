#include <gtest/gtest.h>

#include "../../../include/seqan3/alphabet/quality.hpp"
#include "../../../include/seqan3/alphabet/quality/phred.hpp"

using namespace seqan3;

TEST(illumina18_constructor, default_ctr)
{
    illumina18 illu;
    EXPECT_EQ(illumina18::value_size, 42);
    EXPECT_EQ(illu.offset_char, '!');
}


TEST(illumina18_implicit_assign, implicit_assign)
{
    illumina18 illu;
    illu = 19;
    // expect size unmodified
    EXPECT_EQ(illumina18::value_size, 42);
    // newly assigned member
    EXPECT_EQ(illu.value, 19);
}


TEST(illumina18_op_char, op_char)
{
    illumina18 illu;
    illu = 0;
    char c = char(illu);
    EXPECT_EQ(c, '!');
}


TEST(illumina18_op_tochar, op_tochar)
{
    illumina18 illu;
    illu = 0;
    char c = illu.to_char();
    EXPECT_EQ(c, '!');
    illu = 41;
    c = illu.to_char();
    EXPECT_EQ(c, 'J');
}


TEST(illumina18_from_char, from_char)
{
    illumina18 illu;
    illu = illu.from_char('!');
    EXPECT_EQ(0, illu.value);
}


TEST(illumina18_to_integral, to_integral)
{
    illumina18 illu;
    illu = 0;
    illumina18 illu2 = illu.from_char('!');
    EXPECT_EQ(0, illu2.value);
}


TEST(illumina18_from_integral, from_integral)
{
    illumina18 illu;
    illu = 0;
    illumina18 illu2 = illu.from_integral(1);
    EXPECT_EQ(1, illu2.value);
}

TEST(illumina18_from_phred, from_phred)
{
    illumina18 illu;
    illumina18 illu2 = illu.from_phred(5);
    EXPECT_EQ(5, illu2.value);
}

TEST(illumina18_to_phred, from_phred)
{
    illumina18 illu;
    illumina18 illu2 = illu.from_phred(5);
    EXPECT_EQ(5, illu2.value);
}
