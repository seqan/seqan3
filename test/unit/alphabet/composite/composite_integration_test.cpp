// -----------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2020, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2020, Knut Reinert & MPI für molekulare Genetik
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at: https://github.com/seqan/seqan3/blob/master/LICENSE.md
// -----------------------------------------------------------------------------------------------------

#include <gtest/gtest.h>

#include <seqan3/alphabet/all.hpp>

using seqan3::operator""_aa27;
using seqan3::operator""_dna4;
using seqan3::operator""_rna4;

using qualified_dna_phred42 = seqan3::qualified<seqan3::dna4, seqan3::phred42>;
using qualified_gapped_dna_phred42 = seqan3::qualified<seqan3::gapped<seqan3::dna4>, seqan3::phred42>;
using gapped_qualified_dna_phred42 = seqan3::gapped<qualified_dna_phred42>;
using qualified_qualified_gapped_dna_phred42_phred42 = seqan3::qualified<qualified_gapped_dna_phred42, seqan3::phred42>;
using gapped_alphabet_variant_dna_phred42 = seqan3::gapped<seqan3::alphabet_variant<seqan3::dna4, seqan3::phred42>>;

// Some haessllihckeiten-tests

TEST(composite, custom_constructors)
{
    qualified_dna_phred42 t11{'C'_dna4};
    qualified_dna_phred42 t12{'C'_rna4};
    qualified_dna_phred42 t13{seqan3::phred42{3}};
    qualified_dna_phred42 t14{seqan3::phred63{3}};

    seqan3::qualified<seqan3::aa27, seqan3::phred63> t20{'K'_aa27, seqan3::phred63{}};
    seqan3::qualified<seqan3::aa27, seqan3::phred63> t21{'K'_aa27};
    seqan3::qualified<seqan3::aa27, seqan3::phred63> t22{seqan3::phred63{3}};
    seqan3::qualified<seqan3::aa27, seqan3::phred63> t23{seqan3::phred42{3}};

    qualified_gapped_dna_phred42 t31{'C'_dna4};
    qualified_gapped_dna_phred42 t32{'C'_rna4};
    qualified_gapped_dna_phred42 t33{seqan3::phred42{3}};
    qualified_gapped_dna_phred42 t34{seqan3::gap{}};
    qualified_gapped_dna_phred42 t35{seqan3::gapped<seqan3::dna4>('C'_dna4)};
    qualified_gapped_dna_phred42 t36{seqan3::gapped<seqan3::dna4>(seqan3::gap{})};
    qualified_gapped_dna_phred42 t37{seqan3::gap{}, seqan3::phred42{3}};

    gapped_qualified_dna_phred42 t41{'C'_dna4};
    gapped_qualified_dna_phred42 t42{'C'_rna4};
    gapped_qualified_dna_phred42 t43{seqan3::phred42{3}};
    gapped_qualified_dna_phred42 t44{seqan3::gap{}};
    gapped_qualified_dna_phred42 t45{qualified_dna_phred42{'C'_dna4, seqan3::phred42{0}}};

    qualified_qualified_gapped_dna_phred42_phred42 t51{'C'_dna4};
    qualified_qualified_gapped_dna_phred42_phred42 t52{'C'_rna4};
    qualified_qualified_gapped_dna_phred42_phred42 t53{seqan3::phred42{3}};
    qualified_qualified_gapped_dna_phred42_phred42 t54{seqan3::gap{}};
    qualified_qualified_gapped_dna_phred42_phred42 t55{seqan3::gapped<seqan3::dna4>('C'_dna4)};
    qualified_qualified_gapped_dna_phred42_phred42 t56{seqan3::gapped<seqan3::dna4>(seqan3::gap{})};

    gapped_alphabet_variant_dna_phred42 t61{'C'_dna4};
    gapped_alphabet_variant_dna_phred42 t62{'C'_rna4};
    gapped_alphabet_variant_dna_phred42 t63{seqan3::phred42{3}};
    gapped_alphabet_variant_dna_phred42 t64{seqan3::gap{}};
    gapped_alphabet_variant_dna_phred42 t65{qualified_dna_phred42{'C'_dna4, seqan3::phred42{0}}};

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
    constexpr qualified_dna_phred42 t11{'C'_dna4};
    constexpr qualified_dna_phred42 t12{'C'_rna4};
    constexpr qualified_dna_phred42 t13{seqan3::phred42{3}};
    constexpr qualified_dna_phred42 t14{seqan3::phred63{3}};

    constexpr seqan3::qualified<seqan3::aa27, seqan3::phred63> t21{'K'_aa27};
    constexpr seqan3::qualified<seqan3::aa27, seqan3::phred63> t22{seqan3::phred63{3}};
    constexpr seqan3::qualified<seqan3::aa27, seqan3::phred63> t23{seqan3::phred42{3}};

    constexpr qualified_gapped_dna_phred42 t31{'C'_dna4};
    constexpr qualified_gapped_dna_phred42 t32{'C'_rna4};
    constexpr qualified_gapped_dna_phred42 t33{seqan3::phred42{3}};
    constexpr qualified_gapped_dna_phred42 t34{seqan3::gap{}};
    constexpr qualified_gapped_dna_phred42 t35{seqan3::gapped<seqan3::dna4>('C'_dna4)};
    constexpr qualified_gapped_dna_phred42 t36{seqan3::gapped<seqan3::dna4>(seqan3::gap{})};
    constexpr qualified_gapped_dna_phred42 t37{seqan3::gap{}, seqan3::phred42{3}};

    constexpr gapped_qualified_dna_phred42 t41{'C'_dna4};
    constexpr gapped_qualified_dna_phred42 t42{'C'_rna4};
    constexpr gapped_qualified_dna_phred42 t43{seqan3::phred42{3}};
    constexpr gapped_qualified_dna_phred42 t44{seqan3::gap{}};
    constexpr gapped_qualified_dna_phred42 t45{qualified_dna_phred42{'C'_dna4, seqan3::phred42{0}}};

    constexpr qualified_qualified_gapped_dna_phred42_phred42 t51{'C'_dna4};
    constexpr qualified_qualified_gapped_dna_phred42_phred42 t52{'C'_rna4};
    constexpr qualified_qualified_gapped_dna_phred42_phred42 t53{seqan3::phred42{3}};
    constexpr qualified_qualified_gapped_dna_phred42_phred42 t54{seqan3::gap{}};
    constexpr qualified_qualified_gapped_dna_phred42_phred42 t55{seqan3::gapped<seqan3::dna4>('C'_dna4)};
    constexpr qualified_qualified_gapped_dna_phred42_phred42 t56{seqan3::gapped<seqan3::dna4>(seqan3::gap{})};

    constexpr gapped_alphabet_variant_dna_phred42 t61{'C'_dna4};
    constexpr gapped_alphabet_variant_dna_phred42 t62{'C'_rna4};
    constexpr gapped_alphabet_variant_dna_phred42 t63{seqan3::phred42{3}};
    constexpr gapped_alphabet_variant_dna_phred42 t64{seqan3::gap{}};
    constexpr gapped_alphabet_variant_dna_phred42 t65{qualified_dna_phred42{'C'_dna4, seqan3::phred42{0}}};
}

TEST(composite, custom_assignment)
{
    qualified_dna_phred42 t11{};
    qualified_dna_phred42 t12{'C'_dna4};
    qualified_dna_phred42 t13{'C'_dna4, seqan3::phred42{3}};
    t11 = 'C'_dna4;
    EXPECT_EQ(t11, t12);
    t11 = 'C'_rna4;
    EXPECT_EQ(t11, t12);
    t11 = seqan3::phred42{3};
    EXPECT_EQ(t11, t13);
    // t11 = seqan3::phred63{3}; // does not work because of explicit conversion

    seqan3::qualified<seqan3::aa27, seqan3::phred63> t20{'K'_aa27, seqan3::phred63{}};
    seqan3::qualified<seqan3::aa27, seqan3::phred63> t21{};
    seqan3::qualified<seqan3::aa27, seqan3::phred63> t22{'K'_aa27, seqan3::phred63{3}};
    t21 = 'K'_aa27;
    EXPECT_EQ(t20, t21);
    t21 = seqan3::phred63{3};
    EXPECT_EQ(t21, t22);

    qualified_gapped_dna_phred42 t31{};
    qualified_gapped_dna_phred42 t32{'C'_dna4};
    qualified_gapped_dna_phred42 t33{'C'_dna4, seqan3::phred42{3}};
    qualified_gapped_dna_phred42 t34{seqan3::gap{}, seqan3::phred42{3}};
    t31 = 'C'_dna4;
    EXPECT_EQ(t31, t32);
    t31 = 'C'_rna4;
    EXPECT_EQ(t31, t32);
    t31 = seqan3::phred42{3};
    EXPECT_EQ(t31, t33);
    t31 = seqan3::gap{};
    EXPECT_EQ(t31, t34);
    t31 = seqan3::gapped<seqan3::dna4>('C'_dna4);
    EXPECT_EQ(t31, t33);
    t31 = seqan3::gapped<seqan3::dna4>(seqan3::gap{});
    EXPECT_EQ(t31, t34);

    gapped_qualified_dna_phred42 t41{};
    gapped_qualified_dna_phred42 t42{'C'_dna4};
    gapped_qualified_dna_phred42 t43{qualified_dna_phred42{'C'_dna4, seqan3::phred42{3}}};
    gapped_qualified_dna_phred42 t44{seqan3::gap{}};
    gapped_qualified_dna_phred42 t45{qualified_dna_phred42{'C'_dna4, seqan3::phred42{0}}};
    t41 = 'C'_dna4;
    EXPECT_EQ(t41, t42);
    t41 = 'C'_rna4;
    EXPECT_EQ(t41, t42);
    t41 = seqan3::phred42{3};
    // EXPECT_EQ(t41, t43); should work intuitively but does not because on assignment the seqan3::qualified object is defaulted
    t41 = seqan3::gap{};
    EXPECT_EQ(t41, t44);
    t41 = qualified_dna_phred42{'C'_dna4, seqan3::phred42{0}};
    EXPECT_EQ(t41, t45);

    qualified_qualified_gapped_dna_phred42_phred42 t51{};
    qualified_qualified_gapped_dna_phred42_phred42 t52{'C'_dna4};
    qualified_qualified_gapped_dna_phred42_phred42 t53{qualified_gapped_dna_phred42{'C'_dna4, seqan3::phred42{0}},
                                                       seqan3::phred42{3}};
    qualified_qualified_gapped_dna_phred42_phred42 t54{qualified_gapped_dna_phred42{seqan3::gap{}, seqan3::phred42{0}},
                                                       seqan3::phred42{3}};
    t51 = 'C'_dna4;
    EXPECT_EQ(t51, t52);
    t51 = 'C'_rna4;
    EXPECT_EQ(t51, t52);
    t51 = seqan3::phred42{3};
    EXPECT_EQ(t51, t53);
    t51 = seqan3::gap{};
    EXPECT_EQ(t51, t54);
    t51 = seqan3::gapped<seqan3::dna4>('C'_dna4);
    EXPECT_EQ(t51, t53);
    t51 = seqan3::gapped<seqan3::dna4>(seqan3::gap{});
    EXPECT_EQ(t51, t54);

    gapped_alphabet_variant_dna_phred42 t61{};
    gapped_alphabet_variant_dna_phred42 t62{'C'_dna4};
    gapped_alphabet_variant_dna_phred42 t63{seqan3::phred42{3}};
    gapped_alphabet_variant_dna_phred42 t64{seqan3::gap{}};
    gapped_alphabet_variant_dna_phred42 t65{qualified_dna_phred42{'C'_dna4, seqan3::phred42{0}}};
    t61 = 'C'_dna4;
    EXPECT_EQ(t61, t62);
    t61 = 'C'_rna4;
    EXPECT_EQ(t61, t62);
    t61 = seqan3::phred42{3};
    EXPECT_EQ(t61, t63);
    t61 = seqan3::gap{};
    EXPECT_EQ(t61, t64);
    t61 = qualified_dna_phred42{'C'_dna4, seqan3::phred42{0}};
    EXPECT_EQ(t61, t65);

}

constexpr bool do_assignment()
{
    qualified_dna_phred42 t11{};
    t11 = 'C'_dna4;
    t11 = 'C'_rna4;
    t11 = seqan3::phred42{3};
    // t11 = seqan3::phred63{3}; // does not work because of explicit conversion

    seqan3::qualified<seqan3::aa27, seqan3::phred63> t21{};
    t21 = 'K'_aa27;
    t21 = seqan3::phred63{3};

    qualified_gapped_dna_phred42 t31{};
    t31 = 'C'_dna4;
    t31 = 'C'_rna4;
    t31 = seqan3::phred42{3};
    t31 = seqan3::gap{};
    t31 = seqan3::gapped<seqan3::dna4>('C'_dna4);
    t31 = seqan3::gapped<seqan3::dna4>(seqan3::gap{});

    gapped_qualified_dna_phred42 t41{};
    t41 = 'C'_dna4;
    t41 = 'C'_rna4;
    t41 = seqan3::phred42{3};
    t41 = seqan3::gap{};
    t41 = qualified_dna_phred42{'C'_dna4, seqan3::phred42{0}};

    qualified_qualified_gapped_dna_phred42_phred42 t51{};
    t51 = 'C'_dna4;
    t51 = 'C'_rna4;
    t51 = seqan3::phred42{3};
    t51 = seqan3::gap{};
    t51 = seqan3::gapped<seqan3::dna4>('C'_dna4);
    t51 = seqan3::gapped<seqan3::dna4>(seqan3::gap{});

    gapped_alphabet_variant_dna_phred42 t61{};
    t61 = 'C'_rna4;
    t61 = seqan3::phred42{3};
    t61 = seqan3::gap{};
    t61 = qualified_dna_phred42{'C'_dna4, seqan3::phred42{0}};

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

    qualified_dna_phred42 t11{'C'_dna4, seqan3::phred42{3}};
    EXPECT_EQ(t11, 'C'_dna4);
    EXPECT_EQ(t11, 'C'_rna4);
    EXPECT_EQ(t11, seqan3::phred42{3});
    EXPECT_LT(t11, 'G'_dna4);
    EXPECT_LT(t11, 'G'_rna4);
    EXPECT_LT(t11, seqan3::phred42{4});

    EXPECT_EQ('C'_dna4, t11);
    EXPECT_EQ('C'_rna4, t11);
    EXPECT_EQ(seqan3::phred42{3}, t11);
    EXPECT_LT('A'_dna4, t11);
    EXPECT_LT('A'_rna4, t11);
    EXPECT_LT(seqan3::phred42{2}, t11);

    seqan3::qualified<seqan3::aa27, seqan3::phred63> t21{'K'_aa27, seqan3::phred63{3}};
    EXPECT_EQ(t21, 'K'_aa27);
    EXPECT_EQ(t21, seqan3::phred63{3});
    EXPECT_LT(t21, 'L'_aa27);
    EXPECT_LT(t21, seqan3::phred63{4});

    EXPECT_EQ('K'_aa27, t21);
    EXPECT_EQ(seqan3::phred63{3}, t21);
    EXPECT_LT('C'_aa27, t21);
    EXPECT_LT(seqan3::phred63{2}, t21);

    qualified_gapped_dna_phred42 t31{'C'_dna4, seqan3::phred42{3}};
    EXPECT_EQ(t31, 'C'_dna4);
    EXPECT_EQ(t31, 'C'_rna4);
    EXPECT_EQ(t31, seqan3::phred42{3});
    EXPECT_NE(t31, seqan3::gap{});
    EXPECT_EQ(t31, seqan3::gapped<seqan3::dna4>('C'_dna4));
    EXPECT_LT(t31, 'G'_dna4); // *
    EXPECT_LT(t31, 'G'_rna4); // *
    EXPECT_LT(t31, seqan3::phred42{4});
    EXPECT_LT(t31, seqan3::gap{}); // *
    EXPECT_LT(t31, seqan3::gapped<seqan3::dna4>('G'_dna4));

    EXPECT_EQ('C'_dna4, t31);
    EXPECT_EQ('C'_rna4, t31);
    EXPECT_EQ(seqan3::phred42{3}, t31);
    EXPECT_NE(seqan3::gap{}, t31);
    EXPECT_EQ(seqan3::gapped<seqan3::dna4>('C'_dna4), t31);
    EXPECT_LT('A'_dna4, t31); // *
    EXPECT_LT('A'_rna4, t31); // *
    EXPECT_LT(seqan3::phred42{2}, t31);
    EXPECT_GT(seqan3::gap{}, t31); // *
    EXPECT_LT(seqan3::gapped<seqan3::dna4>('A'_dna4), t31);

    gapped_qualified_dna_phred42 t41{qualified_dna_phred42{'C'_dna4, seqan3::phred42{3}}};
    EXPECT_EQ(t41, 'C'_dna4);
    EXPECT_EQ(t41, 'C'_rna4);
    EXPECT_EQ(t41, seqan3::phred42{3});
    EXPECT_NE(t41, seqan3::gap{});
    EXPECT_EQ(t41, (qualified_dna_phred42{'C'_dna4, seqan3::phred42{3}}));
    EXPECT_EQ(t41, (gapped_qualified_dna_phred42{qualified_dna_phred42{'C'_dna4, seqan3::phred42{3}}}));
//     EXPECT_LT(t41, 'G'_dna4); // not supposed to work
//     EXPECT_LT(t41, 'G'_rna4); // not supposed to work
//     EXPECT_LT(t41, seqan3::phred42{4}); // would never be LT, because seqan3::dna4 part of tuple defaulted to 'A' on RHS
    EXPECT_LT(t41, seqan3::gap{}); // *
    EXPECT_LT(t41, (qualified_dna_phred42{'G'_dna4, seqan3::phred42{2}})); // *
    EXPECT_LT(t41, (gapped_qualified_dna_phred42{qualified_dna_phred42{'G'_dna4, seqan3::phred42{2}}}));

    EXPECT_EQ('C'_dna4, t41);
    EXPECT_EQ('C'_rna4, t41);
    EXPECT_EQ(seqan3::phred42{3}, t41);
    EXPECT_EQ((qualified_dna_phred42{'C'_dna4, seqan3::phred42{3}}), t41);
    EXPECT_NE(seqan3::gap{}, t41);
//     EXPECT_LT('A'_dna4, t41); // not supposed to work
//     EXPECT_LT('A'_rna4, t41); // not supposed to work
//     EXPECT_LT(seqan3::phred42{2}, t41); // not supposed to work
    EXPECT_LT((qualified_dna_phred42{'A'_dna4, seqan3::phred42{2}}), t41); // *
    EXPECT_GT(seqan3::gap{}, t41); // *

    qualified_qualified_gapped_dna_phred42_phred42 t51{qualified_gapped_dna_phred42{'C'_dna4, seqan3::phred42{3}}};
    EXPECT_EQ(t51, 'C'_dna4);
    EXPECT_EQ(t51, 'C'_rna4);
    EXPECT_NE(t51, seqan3::gap{});
    EXPECT_EQ(t51, seqan3::gapped<seqan3::dna4>('C'_dna4));
    EXPECT_EQ(t51, seqan3::phred42{0}); // "outer" phred element
    EXPECT_EQ(t51, (qualified_gapped_dna_phred42{'C'_dna4, seqan3::phred42{3}}));
//     EXPECT_LT(t51, 'G'_dna4); // not supposed to work
//     EXPECT_LT(t51, 'G'_rna4); // not supposed to work
//     EXPECT_LT(t51, seqan3::gap{}); // not supposed to work
//     EXPECT_LT(t51, seqan3::gapped<seqan3::dna4>('G'_dna4)); // not supposed to work
    EXPECT_LT(t51, seqan3::phred42{1});
    EXPECT_LT(t51, (qualified_gapped_dna_phred42{'C'_dna4, seqan3::phred42{4}}));

    EXPECT_EQ('C'_dna4, t51);
    EXPECT_EQ('C'_rna4, t51);
    EXPECT_NE(seqan3::gap{}, t51);
    EXPECT_EQ(seqan3::gapped<seqan3::dna4>('C'_dna4), t51);
    EXPECT_EQ(seqan3::phred42{0}, t51);
    EXPECT_EQ((qualified_gapped_dna_phred42{'C'_dna4, seqan3::phred42{3}}), t51);
//     EXPECT_LT('A'_dna4, t51); // not supposed to work
//     EXPECT_LT('A'_rna4, t51); // not supposed to work
//     EXPECT_GT(seqan3::gap{}, t51); // not supposed to work
//     EXPECT_LT(seqan3::gapped<seqan3::dna4>('A'_dna4), t51); // not supposed to work
    EXPECT_GT(seqan3::phred42{1}, t51);
    EXPECT_GT((qualified_gapped_dna_phred42{'C'_dna4, seqan3::phred42{4}}), t51);

    gapped_alphabet_variant_dna_phred42 t61{'C'_rna4};
    EXPECT_EQ(t61, 'C'_rna4);
    EXPECT_EQ(t61, 'C'_dna4);
    EXPECT_NE(t61, seqan3::gap{});
    EXPECT_NE(t61, seqan3::phred42{0});
    EXPECT_LT(t61, 'G'_rna4); // *
    EXPECT_LT(t61, 'G'_dna4); // *
    EXPECT_LT(t61, seqan3::gap{}); // *
    EXPECT_LT(t61, seqan3::phred42{1}); // *

    EXPECT_EQ('C'_rna4, t61);
    EXPECT_EQ('C'_dna4, t61);
    EXPECT_NE(seqan3::gap{}, t61);
    EXPECT_NE(seqan3::phred42{0}, t61);
    EXPECT_LT('A'_rna4, t61); // *
    EXPECT_LT('A'_dna4, t61); // *
    EXPECT_GT(seqan3::gap{}, t61); // *
    EXPECT_GT(seqan3::phred42{0}, t61); // *
}

TEST(composite, get_)
{
    using seqan3::get;

    qualified_qualified_gapped_dna_phred42_phred42 t51{qualified_gapped_dna_phred42{'C'_dna4, seqan3::phred42{3}}};
    EXPECT_EQ(get<0>(t51),            'C'_dna4);
    EXPECT_EQ(get<0>(get<0>(t51)),    'C'_dna4);

    EXPECT_EQ(get<0>(t51),            'C'_rna4);
    EXPECT_EQ(get<0>(get<0>(t51)),    'C'_rna4);

    EXPECT_NE(get<0>(t51),            seqan3::gap{});
    EXPECT_NE(get<0>(get<0>(t51)),    seqan3::gap{});

    EXPECT_EQ(get<0>(t51),            seqan3::gapped<seqan3::dna4>('C'_dna4));
    EXPECT_EQ(get<0>(get<0>(t51)),    seqan3::gapped<seqan3::dna4>('C'_dna4));

    EXPECT_NE(get<0>(t51),            seqan3::phred42{0});
}
