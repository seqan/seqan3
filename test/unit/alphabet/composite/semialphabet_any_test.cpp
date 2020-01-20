// -----------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2020, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2020, Knut Reinert & MPI für molekulare Genetik
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at: https://github.com/seqan/seqan3/blob/master/LICENSE.md
// -----------------------------------------------------------------------------------------------------

#include <gtest/gtest.h>

#include "../alphabet_test_template.hpp"
#include "../alphabet_constexpr_test_template.hpp"

#include <seqan3/alphabet/aminoacid/all.hpp>
#include <seqan3/alphabet/composite/semialphabet_any.hpp>

using seqan3::operator""_aa10li;
using seqan3::operator""_aa10murphy;

TEST(semialphabet_any_test, initialise_from_alphabet)
{
    seqan3::semialphabet_any<10> letter0{'A'_aa10li};
    seqan3::semialphabet_any<10> letter1{'B'_aa10li};
    seqan3::semialphabet_any<10> letter2{'C'_aa10li};
    seqan3::semialphabet_any<10> letter3{'F'_aa10li};
    seqan3::semialphabet_any<10> letter4{'G'_aa10li};
    seqan3::semialphabet_any<10> letter5{'H'_aa10li};
    seqan3::semialphabet_any<10> letter6{'I'_aa10li};
    seqan3::semialphabet_any<10> letter7{'J'_aa10li};
    seqan3::semialphabet_any<10> letter8{'K'_aa10li};
    seqan3::semialphabet_any<10> letter9{'P'_aa10li};

    EXPECT_EQ(letter0.to_rank(), 0);
    EXPECT_EQ(letter1.to_rank(), 1);
    EXPECT_EQ(letter2.to_rank(), 2);
    EXPECT_EQ(letter3.to_rank(), 3);
    EXPECT_EQ(letter4.to_rank(), 4);
    EXPECT_EQ(letter5.to_rank(), 5);
    EXPECT_EQ(letter6.to_rank(), 6);
    EXPECT_EQ(letter7.to_rank(), 7);
    EXPECT_EQ(letter8.to_rank(), 8);
    EXPECT_EQ(letter9.to_rank(), 9);

    seqan3::aa10murphy new_letter0{'A'_aa10murphy};
    seqan3::aa10murphy new_letter1{'B'_aa10murphy};
    seqan3::aa10murphy new_letter2{'C'_aa10murphy};
    seqan3::aa10murphy new_letter3{'F'_aa10murphy};
    seqan3::aa10murphy new_letter4{'G'_aa10murphy};
    seqan3::aa10murphy new_letter5{'H'_aa10murphy};
    seqan3::aa10murphy new_letter6{'I'_aa10murphy};
    seqan3::aa10murphy new_letter7{'K'_aa10murphy};
    seqan3::aa10murphy new_letter8{'P'_aa10murphy};
    seqan3::aa10murphy new_letter9{'S'_aa10murphy};

    EXPECT_EQ(letter0.to_rank(), new_letter0.to_rank());
    EXPECT_EQ(letter1.to_rank(), new_letter1.to_rank());
    EXPECT_EQ(letter2.to_rank(), new_letter2.to_rank());
    EXPECT_EQ(letter3.to_rank(), new_letter3.to_rank());
    EXPECT_EQ(letter4.to_rank(), new_letter4.to_rank());
    EXPECT_EQ(letter5.to_rank(), new_letter5.to_rank());
    EXPECT_EQ(letter6.to_rank(), new_letter6.to_rank());
    EXPECT_EQ(letter7.to_rank(), new_letter7.to_rank());
    EXPECT_EQ(letter8.to_rank(), new_letter8.to_rank());
    EXPECT_EQ(letter9.to_rank(), new_letter9.to_rank());
}

TEST(semialphabet_any_test, convert_to_alphabet)
{
    seqan3::semialphabet_any<10> letter0{'A'_aa10li};
    seqan3::semialphabet_any<10> letter1{'B'_aa10li};
    seqan3::semialphabet_any<10> letter2{'C'_aa10li};
    seqan3::semialphabet_any<10> letter3{'F'_aa10li};
    seqan3::semialphabet_any<10> letter4{'G'_aa10li};
    seqan3::semialphabet_any<10> letter5{'H'_aa10li};
    seqan3::semialphabet_any<10> letter6{'I'_aa10li};
    seqan3::semialphabet_any<10> letter7{'J'_aa10li};
    seqan3::semialphabet_any<10> letter8{'K'_aa10li};
    seqan3::semialphabet_any<10> letter9{'P'_aa10li};

    EXPECT_EQ(static_cast<seqan3::aa10li>(letter0), 'A'_aa10li);
    EXPECT_EQ(static_cast<seqan3::aa10li>(letter1), 'B'_aa10li);
    EXPECT_EQ(static_cast<seqan3::aa10li>(letter2), 'C'_aa10li);
    EXPECT_EQ(static_cast<seqan3::aa10li>(letter3), 'F'_aa10li);
    EXPECT_EQ(static_cast<seqan3::aa10li>(letter4), 'G'_aa10li);
    EXPECT_EQ(static_cast<seqan3::aa10li>(letter5), 'H'_aa10li);
    EXPECT_EQ(static_cast<seqan3::aa10li>(letter6), 'I'_aa10li);
    EXPECT_EQ(static_cast<seqan3::aa10li>(letter7), 'J'_aa10li);
    EXPECT_EQ(static_cast<seqan3::aa10li>(letter8), 'K'_aa10li);
    EXPECT_EQ(static_cast<seqan3::aa10li>(letter9), 'P'_aa10li);

    EXPECT_EQ(static_cast<seqan3::aa10murphy>(letter0), 'A'_aa10murphy);
    EXPECT_EQ(static_cast<seqan3::aa10murphy>(letter1), 'B'_aa10murphy);
    EXPECT_EQ(static_cast<seqan3::aa10murphy>(letter2), 'C'_aa10murphy);
    EXPECT_EQ(static_cast<seqan3::aa10murphy>(letter3), 'F'_aa10murphy);
    EXPECT_EQ(static_cast<seqan3::aa10murphy>(letter4), 'G'_aa10murphy);
    EXPECT_EQ(static_cast<seqan3::aa10murphy>(letter5), 'H'_aa10murphy);
    EXPECT_EQ(static_cast<seqan3::aa10murphy>(letter6), 'I'_aa10murphy);
    EXPECT_EQ(static_cast<seqan3::aa10murphy>(letter7), 'K'_aa10murphy);
    EXPECT_EQ(static_cast<seqan3::aa10murphy>(letter8), 'P'_aa10murphy);
    EXPECT_EQ(static_cast<seqan3::aa10murphy>(letter9), 'S'_aa10murphy);
}
