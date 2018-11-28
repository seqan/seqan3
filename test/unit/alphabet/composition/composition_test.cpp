// -----------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2018, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2018, Knut Reinert & MPI für molekulare Genetik
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at: https://github.com/seqan/seqan3/blob/master/LICENSE
// -----------------------------------------------------------------------------------------------------

#include <gtest/gtest.h>

#include <seqan3/alphabet/all.hpp>

using namespace seqan3;

// Some haessllihckeiten-tests

TEST(composition, custom_constructors)
{
    qualified<dna4, phred42> t11{dna4::C};
    qualified<dna4, phred42> t12{rna4::C};
    qualified<dna4, phred42> t13{phred42{3}};
    qualified<dna4, phred42> t14{phred63{3}};

    qualified<aa27, phred63> t20{aa27::K, phred63{}};
    qualified<aa27, phred63> t21{aa27::K};
    qualified<aa27, phred63> t22{phred63{3}};
    qualified<aa27, phred63> t23{phred42{3}};

    qualified<gapped<dna4>, phred42> t31{dna4::C};
    qualified<gapped<dna4>, phred42> t32{rna4::C};
    qualified<gapped<dna4>, phred42> t33{phred42{3}};
    qualified<gapped<dna4>, phred42> t34{gap::GAP};
    qualified<gapped<dna4>, phred42> t35{gapped<dna4>(dna4::C)};
    qualified<gapped<dna4>, phred42> t36{gapped<dna4>(gap::GAP)};
    qualified<gapped<dna4>, phred42> t37{gap::GAP, phred42{3}};

    gapped<qualified<dna4, phred42>> t41{dna4::C};
    gapped<qualified<dna4, phred42>> t42{rna4::C};
    gapped<qualified<dna4, phred42>> t43{phred42{3}};
    gapped<qualified<dna4, phred42>> t44{gap::GAP};
    gapped<qualified<dna4, phred42>> t45{qualified<dna4, phred42>{dna4::C, phred42{0}}};

    qualified<qualified<gapped<dna4>, phred42>, phred42> t51{dna4::C};
    qualified<qualified<gapped<dna4>, phred42>, phred42> t52{rna4::C};
    qualified<qualified<gapped<dna4>, phred42>, phred42> t53{phred42{3}};
    qualified<qualified<gapped<dna4>, phred42>, phred42> t54{gap::GAP};
    qualified<qualified<gapped<dna4>, phred42>, phred42> t55{gapped<dna4>(dna4::C)};
    qualified<qualified<gapped<dna4>, phred42>, phred42> t56{gapped<dna4>(gap::GAP)};

    gapped<union_composition<dna4, phred42>> t61{dna4::C};
    gapped<union_composition<dna4, phred42>> t62{rna4::C};
    gapped<union_composition<dna4, phred42>> t63{phred42{3}};
    gapped<union_composition<dna4, phred42>> t64{gap::GAP};
    gapped<union_composition<dna4, phred42>> t65{qualified<dna4, phred42>{dna4::C, phred42{0}}};

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

TEST(composition_constexpr, custom_constructor)
{
    constexpr qualified<dna4, phred42> t11{dna4::C};
    constexpr qualified<dna4, phred42> t12{rna4::C};
    constexpr qualified<dna4, phred42> t13{phred42{3}};
    constexpr qualified<dna4, phred42> t14{phred63{3}};

    constexpr qualified<aa27, phred63> t21{aa27::K};
    constexpr qualified<aa27, phred63> t22{phred63{3}};
    constexpr qualified<aa27, phred63> t23{phred42{3}};

    constexpr qualified<gapped<dna4>, phred42> t31{dna4::C};
    constexpr qualified<gapped<dna4>, phred42> t32{rna4::C};
    constexpr qualified<gapped<dna4>, phred42> t33{phred42{3}};
    constexpr qualified<gapped<dna4>, phred42> t34{gap::GAP};
    constexpr qualified<gapped<dna4>, phred42> t35{gapped<dna4>(dna4::C)};
    constexpr qualified<gapped<dna4>, phred42> t36{gapped<dna4>(gap::GAP)};
    constexpr qualified<gapped<dna4>, phred42> t37{gap::GAP, phred42{3}};

    constexpr gapped<qualified<dna4, phred42>> t41{dna4::C};
    constexpr gapped<qualified<dna4, phred42>> t42{rna4::C};
    constexpr gapped<qualified<dna4, phred42>> t43{phred42{3}};
    constexpr gapped<qualified<dna4, phred42>> t44{gap::GAP};
    constexpr gapped<qualified<dna4, phred42>> t45{qualified<dna4, phred42>{dna4::C, phred42{0}}};

    constexpr qualified<qualified<gapped<dna4>, phred42>, phred42> t51{dna4::C};
    constexpr qualified<qualified<gapped<dna4>, phred42>, phred42> t52{rna4::C};
    constexpr qualified<qualified<gapped<dna4>, phred42>, phred42> t53{phred42{3}};
    constexpr qualified<qualified<gapped<dna4>, phred42>, phred42> t54{gap::GAP};
    constexpr qualified<qualified<gapped<dna4>, phred42>, phred42> t55{gapped<dna4>(dna4::C)};
    constexpr qualified<qualified<gapped<dna4>, phred42>, phred42> t56{gapped<dna4>(gap::GAP)};

    constexpr gapped<union_composition<dna4, phred42>> t61{dna4::C};
    constexpr gapped<union_composition<dna4, phred42>> t62{rna4::C};
    constexpr gapped<union_composition<dna4, phred42>> t63{phred42{3}};
    constexpr gapped<union_composition<dna4, phred42>> t64{gap::GAP};
    constexpr gapped<union_composition<dna4, phred42>> t65{qualified<dna4, phred42>{dna4::C, phred42{0}}};
}

TEST(composition, custom_assignment)
{
    qualified<dna4, phred42> t11{};
    qualified<dna4, phred42> t12{dna4::C};
    qualified<dna4, phred42> t13{dna4::C, phred42{3}};
    t11 = dna4::C;
    EXPECT_EQ(t11, t12);
    t11 = rna4::C;
    EXPECT_EQ(t11, t12);
    t11 = phred42{3};
    EXPECT_EQ(t11, t13);
    // t11 = phred63{3}; // does not work because of explicit conversion

    qualified<aa27, phred63> t20{aa27::K, phred63{}};
    qualified<aa27, phred63> t21{};
    qualified<aa27, phred63> t22{aa27::K, phred63{3}};
    t21 = aa27::K;
    EXPECT_EQ(t20, t21);
    t21 = phred63{3};
    EXPECT_EQ(t21, t22);

    qualified<gapped<dna4>, phred42> t31{};
    qualified<gapped<dna4>, phred42> t32{dna4::C};
    qualified<gapped<dna4>, phred42> t33{dna4::C, phred42{3}};
    qualified<gapped<dna4>, phred42> t34{gap::GAP, phred42{3}};
    t31 = dna4::C;
    EXPECT_EQ(t31, t32);
    t31 = rna4::C;
    EXPECT_EQ(t31, t32);
    t31 = phred42{3};
    EXPECT_EQ(t31, t33);
    t31 = gap::GAP;
    EXPECT_EQ(t31, t34);
    t31 = gapped<dna4>(dna4::C);
    EXPECT_EQ(t31, t33);
    t31 = gapped<dna4>(gap::GAP);
    EXPECT_EQ(t31, t34);

    gapped<qualified<dna4, phred42>> t41{};
    gapped<qualified<dna4, phred42>> t42{dna4::C};
    gapped<qualified<dna4, phred42>> t43{qualified<dna4, phred42>{dna4::C, phred42{3}}};
    gapped<qualified<dna4, phred42>> t44{gap::GAP};
    gapped<qualified<dna4, phred42>> t45{qualified<dna4, phred42>{dna4::C, phred42{0}}};
    t41 = dna4::C;
    EXPECT_EQ(t41, t42);
    t41 = rna4::C;
    EXPECT_EQ(t41, t42);
    t41 = phred42{3};
    // EXPECT_EQ(t41, t43); should work intuitively but does not because on assignment the qualified object is defaulted
    t41 = gap::GAP;
    EXPECT_EQ(t41, t44);
    t41 = qualified<dna4, phred42>{dna4::C, phred42{0}};
    EXPECT_EQ(t41, t45);

    qualified<qualified<gapped<dna4>, phred42>, phred42> t51{};
    qualified<qualified<gapped<dna4>, phred42>, phred42> t52{dna4::C};
    qualified<qualified<gapped<dna4>, phred42>, phred42> t53{qualified<gapped<dna4>, phred42>{dna4::C, phred42{0}}, phred42{3}};
    qualified<qualified<gapped<dna4>, phred42>, phred42> t54{qualified<gapped<dna4>, phred42>{gap::GAP, phred42{0}}, phred42{3}};
    t51 = dna4::C;
    EXPECT_EQ(t51, t52);
    t51 = rna4::C;
    EXPECT_EQ(t51, t52);
    t51 = phred42{3};
    EXPECT_EQ(t51, t53);
    t51 = gap::GAP;
    EXPECT_EQ(t51, t54);
    t51 = gapped<dna4>(dna4::C);
    EXPECT_EQ(t51, t53);
    t51 = gapped<dna4>(gap::GAP);
    EXPECT_EQ(t51, t54);

    gapped<union_composition<dna4, phred42>> t61{};
    gapped<union_composition<dna4, phred42>> t62{dna4::C};
    gapped<union_composition<dna4, phred42>> t63{phred42{3}};
    gapped<union_composition<dna4, phred42>> t64{gap::GAP};
    gapped<union_composition<dna4, phred42>> t65{qualified<dna4, phred42>{dna4::C, phred42{0}}};
    t61 = dna4::C;
    EXPECT_EQ(t61, t62);
    t61 = rna4::C;
    EXPECT_EQ(t61, t62);
    t61 = phred42{3};
    EXPECT_EQ(t61, t63);
    t61 = gap::GAP;
    EXPECT_EQ(t61, t64);
    t61 = qualified<dna4, phred42>{dna4::C, phred42{0}};
    EXPECT_EQ(t61, t65);

}

constexpr bool do_assignment()
{
    qualified<dna4, phred42> t11{};
    t11 = dna4::C;
    t11 = rna4::C;
    t11 = phred42{3};
    // t11 = phred63{3}; // does not work because of explicit conversion

    qualified<aa27, phred63> t21{};
    t21 = aa27::K;
    t21 = phred63{3};

    qualified<gapped<dna4>, phred42> t31{};
    t31 = dna4::C;
    t31 = rna4::C;
    t31 = phred42{3};
    t31 = gap::GAP;
    t31 = gapped<dna4>(dna4::C);
    t31 = gapped<dna4>(gap::GAP);

    gapped<qualified<dna4, phred42>> t41{};
    t41 = dna4::C;
    t41 = rna4::C;
    t41 = phred42{3};
    t41 = gap::GAP;
    t41 = qualified<dna4, phred42>{dna4::C, phred42{0}};

    qualified<qualified<gapped<dna4>, phred42>, phred42> t51{};
    t51 = dna4::C;
    t51 = rna4::C;
    t51 = phred42{3};
    t51 = gap::GAP;
    t51 = gapped<dna4>(dna4::C);
    t51 = gapped<dna4>(gap::GAP);

    gapped<union_composition<dna4, phred42>> t61{};
    t61 = rna4::C;
    t61 = phred42{3};
    t61 = gap::GAP;
    t61 = qualified<dna4, phred42>{dna4::C, phred42{0}};

    return true;
}

TEST(composition_constexpr, custom_assignment)
{
    [[maybe_unused]] constexpr bool foo = do_assignment();
}

TEST(composition, custom_comparison)
{
    qualified<dna4, phred42> t11{dna4::C, phred42{3}};
    EXPECT_EQ(t11, dna4::C);
    EXPECT_EQ(t11, rna4::C);
    EXPECT_EQ(t11, phred42{3});

    EXPECT_EQ(dna4::C,    t11);
    EXPECT_EQ(rna4::C,    t11);
    EXPECT_EQ(phred42{3}, t11);

    qualified<aa27, phred63> t21{aa27::K, phred63{3}};
    EXPECT_EQ(t21, aa27::K);
    EXPECT_EQ(t21, phred63{3});
    EXPECT_EQ(aa27::K,     t21);
    EXPECT_EQ(phred63{3},  t21);

    qualified<gapped<dna4>, phred42> t31{dna4::C, phred42{3}};
    EXPECT_EQ(t31, dna4::C);
    EXPECT_EQ(t31, rna4::C);
    EXPECT_EQ(t31, phred42{3});
    EXPECT_NE(t31, gap::GAP);
    EXPECT_EQ(t31, gapped<dna4>(dna4::C));

    EXPECT_EQ(dna4::C,                t31);
    EXPECT_EQ(rna4::C,                t31);
    EXPECT_EQ(phred42{3},             t31);
    EXPECT_NE(gap::GAP,               t31);
    EXPECT_EQ(gapped<dna4>(dna4::C),  t31);

    gapped<qualified<dna4, phred42>> t41{qualified<dna4, phred42>{dna4::C, phred42{3}}};
    gapped<qualified<dna4, phred42>> t42{qualified<dna4, phred42>{dna4::C, phred42{0}}};

    EXPECT_EQ(t41, (gapped<qualified<dna4, phred42>>{qualified<dna4, phred42>{dna4::C, phred42{3}}}));
    EXPECT_EQ(t42, dna4::C);
    EXPECT_NE(t41, gap::GAP);
    EXPECT_NE(gap::GAP, t41);

    qualified<qualified<gapped<dna4>, phred42>, phred42> t51{qualified<gapped<dna4>, phred42>{dna4::C, phred42{3}}};
    EXPECT_EQ(t51, dna4::C);
    EXPECT_EQ(t51, rna4::C);
    EXPECT_NE(t51, gap::GAP);
    EXPECT_EQ(t51, gapped<dna4>(dna4::C));
    EXPECT_EQ(t51, phred42{0});

    EXPECT_EQ(dna4::C,                t51);
    EXPECT_EQ(rna4::C,                t51);
    EXPECT_EQ(phred42{0},             t51);
    EXPECT_NE(gap::GAP,               t51);
    EXPECT_EQ(gapped<dna4>(dna4::C),  t51);

    gapped<union_composition<dna4, phred42>> t61{rna4::C};
    EXPECT_EQ(t61, rna4::C);
    EXPECT_EQ(t61, dna4::C);
    EXPECT_NE(t61, gap::GAP);
    EXPECT_NE(t61, phred42{0});

    EXPECT_EQ(rna4::C,    t61);
    EXPECT_EQ(dna4::C,    t61);
    EXPECT_NE(gap::GAP,   t61);
    EXPECT_NE(phred42{0}, t61);

    EXPECT_EQ(t41, dna4::C);
    EXPECT_EQ(t41, rna4::C);
    EXPECT_EQ(t41, phred42{3});
    EXPECT_EQ(t41, (qualified<dna4, phred42>{dna4::C, phred42{3}}));

    EXPECT_EQ(dna4::C,                                         t41);
    EXPECT_EQ(rna4::C,                                         t41);
    EXPECT_EQ(phred42{3},                                      t41);
    EXPECT_EQ((qualified<dna4, phred42>{dna4::C, phred42{3}}), t41);

    EXPECT_EQ(t51, (qualified<gapped<dna4>, phred42>{dna4::C, phred42{3}}));
    EXPECT_EQ((qualified<dna4, phred42>{dna4::C, phred42{3}}), t51);
}
