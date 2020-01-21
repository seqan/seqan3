// -----------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2020, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2020, Knut Reinert & MPI für molekulare Genetik
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at: https://github.com/seqan/seqan3/blob/master/LICENSE.md
// -----------------------------------------------------------------------------------------------------

#include <gtest/gtest.h>

#include <seqan3/alphabet/all.hpp>

using namespace seqan3;

// Some haessllihckeiten-tests

TEST(composite, custom_constructors)
{
    qualified<dna4, phred42> t11{'C'_dna4};
    qualified<dna4, phred42> t12{'C'_rna4};
    qualified<dna4, phred42> t13{phred42{3}};
    qualified<dna4, phred42> t14{phred63{3}};

    qualified<aa27, phred63> t20{'K'_aa27, phred63{}};
    qualified<aa27, phred63> t21{'K'_aa27};
    qualified<aa27, phred63> t22{phred63{3}};
    qualified<aa27, phred63> t23{phred42{3}};

    qualified<gapped<dna4>, phred42> t31{'C'_dna4};
    qualified<gapped<dna4>, phred42> t32{'C'_rna4};
    qualified<gapped<dna4>, phred42> t33{phred42{3}};
    qualified<gapped<dna4>, phred42> t34{gap{}};
    qualified<gapped<dna4>, phred42> t35{gapped<dna4>('C'_dna4)};
    qualified<gapped<dna4>, phred42> t36{gapped<dna4>(gap{})};
    qualified<gapped<dna4>, phred42> t37{gap{}, phred42{3}};

    gapped<qualified<dna4, phred42>> t41{'C'_dna4};
    gapped<qualified<dna4, phred42>> t42{'C'_rna4};
    gapped<qualified<dna4, phred42>> t43{phred42{3}};
    gapped<qualified<dna4, phred42>> t44{gap{}};
    gapped<qualified<dna4, phred42>> t45{qualified<dna4, phred42>{'C'_dna4, phred42{0}}};

    qualified<qualified<gapped<dna4>, phred42>, phred42> t51{'C'_dna4};
    qualified<qualified<gapped<dna4>, phred42>, phred42> t52{'C'_rna4};
    qualified<qualified<gapped<dna4>, phred42>, phred42> t53{phred42{3}};
    qualified<qualified<gapped<dna4>, phred42>, phred42> t54{gap{}};
    qualified<qualified<gapped<dna4>, phred42>, phred42> t55{gapped<dna4>('C'_dna4)};
    qualified<qualified<gapped<dna4>, phred42>, phred42> t56{gapped<dna4>(gap{})};

    gapped<alphabet_variant<dna4, phred42>> t61{'C'_dna4};
    gapped<alphabet_variant<dna4, phred42>> t62{'C'_rna4};
    gapped<alphabet_variant<dna4, phred42>> t63{phred42{3}};
    gapped<alphabet_variant<dna4, phred42>> t64{gap{}};
    gapped<alphabet_variant<dna4, phred42>> t65{qualified<dna4, phred42>{'C'_dna4, phred42{0}}};

    EXPECT_EQ(t11, t12);
    EXPECT_EQ(t13, t14);

    EXPECT_EQ(t20, t21);
    EXPECT_EQ(t22, t23);

    EXPECT_EQ(t31, t32);
    EXPECT_NE(t31, t33);
    EXPECT_NE(t31, t34);
    EXPECT_EQ(t31, t35);
    EXPECT_EQ(t34, t36);

    EXPECT_EQ(t41, t42);
    EXPECT_NE(t41, t43);
    EXPECT_NE(t41, t44);
    EXPECT_EQ(t41, t45);

    EXPECT_EQ(t51, t52);
    EXPECT_NE(t51, t53);
    EXPECT_NE(t51, t54);
    EXPECT_EQ(t51, t55);
    EXPECT_EQ(t54, t56);

    EXPECT_EQ(t61, t62);
    EXPECT_NE(t61, t63);
    EXPECT_NE(t61, t64);
    EXPECT_EQ(t61, t65);
}

TEST(composite_constexpr, custom_constructor)
{
    constexpr qualified<dna4, phred42> t11{'C'_dna4};
    constexpr qualified<dna4, phred42> t12{'C'_rna4};
    constexpr qualified<dna4, phred42> t13{phred42{3}};
    constexpr qualified<dna4, phred42> t14{phred63{3}};

    constexpr qualified<aa27, phred63> t21{'K'_aa27};
    constexpr qualified<aa27, phred63> t22{phred63{3}};
    constexpr qualified<aa27, phred63> t23{phred42{3}};

    constexpr qualified<gapped<dna4>, phred42> t31{'C'_dna4};
    constexpr qualified<gapped<dna4>, phred42> t32{'C'_rna4};
    constexpr qualified<gapped<dna4>, phred42> t33{phred42{3}};
    constexpr qualified<gapped<dna4>, phred42> t34{gap{}};
    constexpr qualified<gapped<dna4>, phred42> t35{gapped<dna4>('C'_dna4)};
    constexpr qualified<gapped<dna4>, phred42> t36{gapped<dna4>(gap{})};
    constexpr qualified<gapped<dna4>, phred42> t37{gap{}, phred42{3}};

    constexpr gapped<qualified<dna4, phred42>> t41{'C'_dna4};
    constexpr gapped<qualified<dna4, phred42>> t42{'C'_rna4};
    constexpr gapped<qualified<dna4, phred42>> t43{phred42{3}};
    constexpr gapped<qualified<dna4, phred42>> t44{gap{}};
    constexpr gapped<qualified<dna4, phred42>> t45{qualified<dna4, phred42>{'C'_dna4, phred42{0}}};

    constexpr qualified<qualified<gapped<dna4>, phred42>, phred42> t51{'C'_dna4};
    constexpr qualified<qualified<gapped<dna4>, phred42>, phred42> t52{'C'_rna4};
    constexpr qualified<qualified<gapped<dna4>, phred42>, phred42> t53{phred42{3}};
    constexpr qualified<qualified<gapped<dna4>, phred42>, phred42> t54{gap{}};
    constexpr qualified<qualified<gapped<dna4>, phred42>, phred42> t55{gapped<dna4>('C'_dna4)};
    constexpr qualified<qualified<gapped<dna4>, phred42>, phred42> t56{gapped<dna4>(gap{})};

    constexpr gapped<alphabet_variant<dna4, phred42>> t61{'C'_dna4};
    constexpr gapped<alphabet_variant<dna4, phred42>> t62{'C'_rna4};
    constexpr gapped<alphabet_variant<dna4, phred42>> t63{phred42{3}};
    constexpr gapped<alphabet_variant<dna4, phred42>> t64{gap{}};
    constexpr gapped<alphabet_variant<dna4, phred42>> t65{qualified<dna4, phred42>{'C'_dna4, phred42{0}}};
}

TEST(composite, custom_assignment)
{
    qualified<dna4, phred42> t11{};
    qualified<dna4, phred42> t12{'C'_dna4};
    qualified<dna4, phred42> t13{'C'_dna4, phred42{3}};
    t11 = 'C'_dna4;
    EXPECT_EQ(t11, t12);
    t11 = 'C'_rna4;
    EXPECT_EQ(t11, t12);
    t11 = phred42{3};
    EXPECT_EQ(t11, t13);
    // t11 = phred63{3}; // does not work because of explicit conversion

    qualified<aa27, phred63> t20{'K'_aa27, phred63{}};
    qualified<aa27, phred63> t21{};
    qualified<aa27, phred63> t22{'K'_aa27, phred63{3}};
    t21 = 'K'_aa27;
    EXPECT_EQ(t20, t21);
    t21 = phred63{3};
    EXPECT_EQ(t21, t22);

    qualified<gapped<dna4>, phred42> t31{};
    qualified<gapped<dna4>, phred42> t32{'C'_dna4};
    qualified<gapped<dna4>, phred42> t33{'C'_dna4, phred42{3}};
    qualified<gapped<dna4>, phred42> t34{gap{}, phred42{3}};
    t31 = 'C'_dna4;
    EXPECT_EQ(t31, t32);
    t31 = 'C'_rna4;
    EXPECT_EQ(t31, t32);
    t31 = phred42{3};
    EXPECT_EQ(t31, t33);
    t31 = gap{};
    EXPECT_EQ(t31, t34);
    t31 = gapped<dna4>('C'_dna4);
    EXPECT_EQ(t31, t33);
    t31 = gapped<dna4>(gap{});
    EXPECT_EQ(t31, t34);

    gapped<qualified<dna4, phred42>> t41{};
    gapped<qualified<dna4, phred42>> t42{'C'_dna4};
    gapped<qualified<dna4, phred42>> t43{qualified<dna4, phred42>{'C'_dna4, phred42{3}}};
    gapped<qualified<dna4, phred42>> t44{gap{}};
    gapped<qualified<dna4, phred42>> t45{qualified<dna4, phred42>{'C'_dna4, phred42{0}}};
    t41 = 'C'_dna4;
    EXPECT_EQ(t41, t42);
    t41 = 'C'_rna4;
    EXPECT_EQ(t41, t42);
    t41 = phred42{3};
    // EXPECT_EQ(t41, t43); should work intuitively but does not because on assignment the qualified object is defaulted
    t41 = gap{};
    EXPECT_EQ(t41, t44);
    t41 = qualified<dna4, phred42>{'C'_dna4, phred42{0}};
    EXPECT_EQ(t41, t45);

    qualified<qualified<gapped<dna4>, phred42>, phred42> t51{};
    qualified<qualified<gapped<dna4>, phred42>, phred42> t52{'C'_dna4};
    qualified<qualified<gapped<dna4>, phred42>, phred42> t53{qualified<gapped<dna4>, phred42>{'C'_dna4, phred42{0}}, phred42{3}};
    qualified<qualified<gapped<dna4>, phred42>, phred42> t54{qualified<gapped<dna4>, phred42>{gap{}, phred42{0}}, phred42{3}};
    t51 = 'C'_dna4;
    EXPECT_EQ(t51, t52);
    t51 = 'C'_rna4;
    EXPECT_EQ(t51, t52);
    t51 = phred42{3};
    EXPECT_EQ(t51, t53);
    t51 = gap{};
    EXPECT_EQ(t51, t54);
    t51 = gapped<dna4>('C'_dna4);
    EXPECT_EQ(t51, t53);
    t51 = gapped<dna4>(gap{});
    EXPECT_EQ(t51, t54);

    gapped<alphabet_variant<dna4, phred42>> t61{};
    gapped<alphabet_variant<dna4, phred42>> t62{'C'_dna4};
    gapped<alphabet_variant<dna4, phred42>> t63{phred42{3}};
    gapped<alphabet_variant<dna4, phred42>> t64{gap{}};
    gapped<alphabet_variant<dna4, phred42>> t65{qualified<dna4, phred42>{'C'_dna4, phred42{0}}};
    t61 = 'C'_dna4;
    EXPECT_EQ(t61, t62);
    t61 = 'C'_rna4;
    EXPECT_EQ(t61, t62);
    t61 = phred42{3};
    EXPECT_EQ(t61, t63);
    t61 = gap{};
    EXPECT_EQ(t61, t64);
    t61 = qualified<dna4, phred42>{'C'_dna4, phred42{0}};
    EXPECT_EQ(t61, t65);

}

constexpr bool do_assignment()
{
    qualified<dna4, phred42> t11{};
    t11 = 'C'_dna4;
    t11 = 'C'_rna4;
    t11 = phred42{3};
    // t11 = phred63{3}; // does not work because of explicit conversion

    qualified<aa27, phred63> t21{};
    t21 = 'K'_aa27;
    t21 = phred63{3};

    qualified<gapped<dna4>, phred42> t31{};
    t31 = 'C'_dna4;
    t31 = 'C'_rna4;
    t31 = phred42{3};
    t31 = gap{};
    t31 = gapped<dna4>('C'_dna4);
    t31 = gapped<dna4>(gap{});

    gapped<qualified<dna4, phred42>> t41{};
    t41 = 'C'_dna4;
    t41 = 'C'_rna4;
    t41 = phred42{3};
    t41 = gap{};
    t41 = qualified<dna4, phred42>{'C'_dna4, phred42{0}};

    qualified<qualified<gapped<dna4>, phred42>, phred42> t51{};
    t51 = 'C'_dna4;
    t51 = 'C'_rna4;
    t51 = phred42{3};
    t51 = gap{};
    t51 = gapped<dna4>('C'_dna4);
    t51 = gapped<dna4>(gap{});

    gapped<alphabet_variant<dna4, phred42>> t61{};
    t61 = 'C'_rna4;
    t61 = phred42{3};
    t61 = gap{};
    t61 = qualified<dna4, phred42>{'C'_dna4, phred42{0}};

    return true;
}

TEST(composite_constexpr, custom_assignment)
{
    [[maybe_unused]] constexpr bool foo = do_assignment();
}

TEST(composite, custom_comparison)
{
    /* Tests marked with "// *" would not be possible if all single argument constructors of alphabet_variant
     * are made explicit */

    qualified<dna4, phred42> t11{'C'_dna4, phred42{3}};
    EXPECT_EQ(t11, 'C'_dna4);
    EXPECT_EQ(t11, 'C'_rna4);
    EXPECT_EQ(t11, phred42{3});
    EXPECT_LT(t11, 'G'_dna4);
    EXPECT_LT(t11, 'G'_rna4);
    EXPECT_LT(t11, phred42{4});

    EXPECT_EQ('C'_dna4,    t11);
    EXPECT_EQ('C'_rna4,    t11);
    EXPECT_EQ(phred42{3},  t11);
    EXPECT_LT('A'_dna4,    t11);
    EXPECT_LT('A'_rna4,    t11);
    EXPECT_LT(phred42{2},  t11);

    qualified<aa27, phred63> t21{'K'_aa27, phred63{3}};
    EXPECT_EQ(t21, 'K'_aa27);
    EXPECT_EQ(t21, phred63{3});
    EXPECT_LT(t21, 'L'_aa27);
    EXPECT_LT(t21, phred63{4});

    EXPECT_EQ('K'_aa27,    t21);
    EXPECT_EQ(phred63{3},  t21);
    EXPECT_LT('C'_aa27,    t21);
    EXPECT_LT(phred63{2},  t21);

    qualified<gapped<dna4>, phred42> t31{'C'_dna4, phred42{3}};
    EXPECT_EQ(t31, 'C'_dna4);
    EXPECT_EQ(t31, 'C'_rna4);
    EXPECT_EQ(t31, phred42{3});
    EXPECT_NE(t31, gap{});
    EXPECT_EQ(t31, gapped<dna4>('C'_dna4));
    EXPECT_LT(t31, 'G'_dna4);                     // *
    EXPECT_LT(t31, 'G'_rna4);                     // *
    EXPECT_LT(t31, phred42{4});
    EXPECT_LT(t31, gap{});                        // *
    EXPECT_LT(t31, gapped<dna4>('G'_dna4));

    EXPECT_EQ('C'_dna4,                t31);
    EXPECT_EQ('C'_rna4,                t31);
    EXPECT_EQ(phred42{3},              t31);
    EXPECT_NE(gap{},                   t31);
    EXPECT_EQ(gapped<dna4>('C'_dna4),  t31);
    EXPECT_LT('A'_dna4,                t31);      // *
    EXPECT_LT('A'_rna4,                t31);      // *
    EXPECT_LT(phred42{2},              t31);
    EXPECT_GT(gap{},                   t31);      // *
    EXPECT_LT(gapped<dna4>('A'_dna4),  t31);

    gapped<qualified<dna4, phred42>> t41{qualified<dna4, phred42>{'C'_dna4, phred42{3}}};
    EXPECT_EQ(t41, 'C'_dna4);
    EXPECT_EQ(t41, 'C'_rna4);
    EXPECT_EQ(t41, phred42{3});
    EXPECT_NE(t41, gap{});
    EXPECT_EQ(t41, (qualified<dna4, phred42>{'C'_dna4, phred42{3}}));
    EXPECT_EQ(t41, (gapped<qualified<dna4, phred42>>{qualified<dna4, phred42>{'C'_dna4, phred42{3}}}));
//     EXPECT_LT(t41, 'G'_dna4);       // not supposed to work
//     EXPECT_LT(t41, 'G'_rna4);       // not supposed to work
//     EXPECT_LT(t41, phred42{4});     // would never be LT, because dna4 part of tuple defaulted to 'A' on RHS
    EXPECT_LT(t41, gap{});                                                                                   // *
    EXPECT_LT(t41, (qualified<dna4, phred42>{'G'_dna4, phred42{2}}));                                        // *
    EXPECT_LT(t41, (gapped<qualified<dna4, phred42>>{qualified<dna4, phred42>{'G'_dna4, phred42{2}}}));

    EXPECT_EQ('C'_dna4,                                         t41);
    EXPECT_EQ('C'_rna4,                                         t41);
    EXPECT_EQ(phred42{3},                                       t41);
    EXPECT_EQ((qualified<dna4, phred42>{'C'_dna4, phred42{3}}), t41);
    EXPECT_NE(gap{},                                            t41);
//     EXPECT_LT('A'_dna4,                                         t41);  // not supposed to work
//     EXPECT_LT('A'_rna4,                                         t41);  // not supposed to work
//     EXPECT_LT(phred42{2},                                       t41);  // not supposed to work
    EXPECT_LT((qualified<dna4, phred42>{'A'_dna4, phred42{2}}), t41);  // *
    EXPECT_GT(gap{},                                            t41);  // *

    qualified<qualified<gapped<dna4>, phred42>, phred42> t51{qualified<gapped<dna4>, phred42>{'C'_dna4, phred42{3}}};
    EXPECT_EQ(t51, 'C'_dna4);
    EXPECT_EQ(t51, 'C'_rna4);
    EXPECT_NE(t51, gap{});
    EXPECT_EQ(t51, gapped<dna4>('C'_dna4));
    EXPECT_EQ(t51, phred42{0}); // "outer" phred element
    EXPECT_EQ(t51, (qualified<gapped<dna4>, phred42>{'C'_dna4, phred42{3}}));
//     EXPECT_LT(t51, 'G'_dna4);                                                      // not supposed to work
//     EXPECT_LT(t51, 'G'_rna4);                                                      // not supposed to work
//     EXPECT_LT(t51, gap{});                                                         // not supposed to work
//     EXPECT_LT(t51, gapped<dna4>('G'_dna4));                                        // not supposed to work
    EXPECT_LT(t51, phred42{1});
    EXPECT_LT(t51, (qualified<gapped<dna4>, phred42>{'C'_dna4, phred42{4}}));

    EXPECT_EQ('C'_dna4,                                                 t51);
    EXPECT_EQ('C'_rna4,                                                 t51);
    EXPECT_NE(gap{},                                                    t51);
    EXPECT_EQ(gapped<dna4>('C'_dna4),                                   t51);
    EXPECT_EQ(phred42{0},                                               t51);
    EXPECT_EQ((qualified<gapped<dna4>, phred42>{'C'_dna4, phred42{3}}), t51);
//     EXPECT_LT('A'_dna4,                                                 t51);      // not supposed to work
//     EXPECT_LT('A'_rna4,                                                 t51);      // not supposed to work
//     EXPECT_GT(gap{},                                                    t51);      // not supposed to work
//     EXPECT_LT(gapped<dna4>('A'_dna4),                                   t51);      // not supposed to work
    EXPECT_GT(phred42{1},                                               t51);
    EXPECT_GT((qualified<gapped<dna4>, phred42>{'C'_dna4, phred42{4}}), t51);

    gapped<alphabet_variant<dna4, phred42>> t61{'C'_rna4};
    EXPECT_EQ(t61, 'C'_rna4);
    EXPECT_EQ(t61, 'C'_dna4);
    EXPECT_NE(t61, gap{});
    EXPECT_NE(t61, phred42{0});
    EXPECT_LT(t61, 'G'_rna4);              // *
    EXPECT_LT(t61, 'G'_dna4);              // *
    EXPECT_LT(t61, gap{});                 // *
    EXPECT_LT(t61, phred42{1});            // *

    EXPECT_EQ('C'_rna4,    t61);
    EXPECT_EQ('C'_dna4,    t61);
    EXPECT_NE(gap{},       t61);
    EXPECT_NE(phred42{0},  t61);
    EXPECT_LT('A'_rna4,    t61);           // *
    EXPECT_LT('A'_dna4,    t61);           // *
    EXPECT_GT(gap{},       t61);           // *
    EXPECT_GT(phred42{0},  t61);           // *
}

TEST(composite, get_)
{
    qualified<qualified<gapped<dna4>, phred42>, phred42> t51{qualified<gapped<dna4>, phred42>{'C'_dna4, phred42{3}}};
    EXPECT_EQ(get<0>(t51),            'C'_dna4);
    EXPECT_EQ(get<0>(get<0>(t51)),    'C'_dna4);

    EXPECT_EQ(get<0>(t51),            'C'_rna4);
    EXPECT_EQ(get<0>(get<0>(t51)),    'C'_rna4);

    EXPECT_NE(get<0>(t51),            gap{});
    EXPECT_NE(get<0>(get<0>(t51)),    gap{});

    EXPECT_EQ(get<0>(t51),            gapped<dna4>('C'_dna4));
    EXPECT_EQ(get<0>(get<0>(t51)),    gapped<dna4>('C'_dna4));

    EXPECT_NE(get<0>(t51),            phred42{0});
}
