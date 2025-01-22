// SPDX-FileCopyrightText: 2006-2025 Knut Reinert & Freie Universität Berlin
// SPDX-FileCopyrightText: 2016-2025 Knut Reinert & MPI für molekulare Genetik
// SPDX-License-Identifier: BSD-3-Clause

#include <gtest/gtest.h>

#include <seqan3/alphabet/aminoacid/aa27.hpp>
#include <seqan3/alphabet/gap/gapped.hpp>
#include <seqan3/alphabet/nucleotide/dna4.hpp>
#include <seqan3/alphabet/nucleotide/rna4.hpp>
#include <seqan3/alphabet/quality/phred42.hpp>
#include <seqan3/alphabet/quality/phred63.hpp>
#include <seqan3/alphabet/quality/qualified.hpp>

#include "../alphabet_constexpr_test_template.hpp"
#include "../alphabet_test_template.hpp"
#include "../semi_alphabet_constexpr_test_template.hpp"
#include "../semi_alphabet_test_template.hpp"

using seqan3::operator""_aa27;
using seqan3::operator""_dna4;
using seqan3::operator""_phred42;
using seqan3::operator""_phred63;
using seqan3::operator""_rna4;

// tests various combinations of alphabet_variant and alphabet_tuple
using qualified_dna_phred42 = seqan3::qualified<seqan3::dna4, seqan3::phred42>;
using qualified_gapped_dna_phred42 = seqan3::qualified<seqan3::gapped<seqan3::dna4>, seqan3::phred42>;
using gapped_qualified_dna_phred42 = seqan3::gapped<qualified_dna_phred42>;
using qualified_qualified_gapped_dna_phred42_phred42 = seqan3::qualified<qualified_gapped_dna_phred42, seqan3::phred42>;
using gapped_alphabet_variant_dna_phred42 = seqan3::gapped<seqan3::alphabet_variant<seqan3::dna4, seqan3::phred42>>;

using alphabet_types = ::testing::Types<qualified_dna_phred42,
                                        qualified_gapped_dna_phred42,
                                        gapped_qualified_dna_phred42,
                                        qualified_qualified_gapped_dna_phred42_phred42,
                                        gapped_alphabet_variant_dna_phred42>;

INSTANTIATE_TYPED_TEST_SUITE_P(alphabet_composite_integration, alphabet, alphabet_types, );
INSTANTIATE_TYPED_TEST_SUITE_P(alphabet_composite_integration, semi_alphabet_test, alphabet_types, );
INSTANTIATE_TYPED_TEST_SUITE_P(alphabet_composite_integration, alphabet_constexpr, alphabet_types, );
INSTANTIATE_TYPED_TEST_SUITE_P(alphabet_composite_integration, semi_alphabet_constexpr, alphabet_types, );

TEST(composite, custom_constructors)
{
    qualified_dna_phred42 t11{'C'_dna4};
    qualified_dna_phred42 t12{'C'_rna4};
    qualified_dna_phred42 t13{'$'_phred42};
    qualified_dna_phred42 t14{'$'_phred63};

    seqan3::qualified<seqan3::aa27, seqan3::phred63> t20{'K'_aa27, seqan3::phred63{}};
    seqan3::qualified<seqan3::aa27, seqan3::phred63> t21{'K'_aa27};
    seqan3::qualified<seqan3::aa27, seqan3::phred63> t22{'$'_phred63};
    seqan3::qualified<seqan3::aa27, seqan3::phred63> t23{'$'_phred42};

    qualified_gapped_dna_phred42 t31{'C'_dna4};
    qualified_gapped_dna_phred42 t32{'C'_rna4};
    qualified_gapped_dna_phred42 t33{'$'_phred42};
    qualified_gapped_dna_phred42 t34{seqan3::gap{}};
    qualified_gapped_dna_phred42 t35{seqan3::gapped<seqan3::dna4>('C'_dna4)};
    qualified_gapped_dna_phred42 t36{seqan3::gapped<seqan3::dna4>(seqan3::gap{})};
    [[maybe_unused]] qualified_gapped_dna_phred42 t37{seqan3::gap{}, '$'_phred42};

    gapped_qualified_dna_phred42 t41{'C'_dna4};
    gapped_qualified_dna_phred42 t42{'C'_rna4};
    gapped_qualified_dna_phred42 t43{'$'_phred42};
    gapped_qualified_dna_phred42 t44{seqan3::gap{}};
    gapped_qualified_dna_phred42 t45{qualified_dna_phred42{'C'_dna4, '!'_phred42}};

    qualified_qualified_gapped_dna_phred42_phred42 t51{'C'_dna4};
    qualified_qualified_gapped_dna_phred42_phred42 t52{'C'_rna4};
    qualified_qualified_gapped_dna_phred42_phred42 t53{'$'_phred42};
    qualified_qualified_gapped_dna_phred42_phred42 t54{seqan3::gap{}};
    qualified_qualified_gapped_dna_phred42_phred42 t55{seqan3::gapped<seqan3::dna4>('C'_dna4)};
    qualified_qualified_gapped_dna_phred42_phred42 t56{seqan3::gapped<seqan3::dna4>(seqan3::gap{})};

    gapped_alphabet_variant_dna_phred42 t61{'C'_dna4};
    gapped_alphabet_variant_dna_phred42 t62{'C'_rna4};
    gapped_alphabet_variant_dna_phred42 t63{'$'_phred42};
    gapped_alphabet_variant_dna_phred42 t64{seqan3::gap{}};
    gapped_alphabet_variant_dna_phred42 t65{qualified_dna_phred42{'C'_dna4, '!'_phred42}};

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
    [[maybe_unused]] constexpr qualified_dna_phred42 t11{'C'_dna4};
    [[maybe_unused]] constexpr qualified_dna_phred42 t12{'C'_rna4};
    [[maybe_unused]] constexpr qualified_dna_phred42 t13{'$'_phred42};
    [[maybe_unused]] constexpr qualified_dna_phred42 t14{'$'_phred63};

    [[maybe_unused]] constexpr seqan3::qualified<seqan3::aa27, seqan3::phred63> t21{'K'_aa27};
    [[maybe_unused]] constexpr seqan3::qualified<seqan3::aa27, seqan3::phred63> t22{'$'_phred63};
    [[maybe_unused]] constexpr seqan3::qualified<seqan3::aa27, seqan3::phred63> t23{'$'_phred42};

    [[maybe_unused]] constexpr qualified_gapped_dna_phred42 t31{'C'_dna4};
    [[maybe_unused]] constexpr qualified_gapped_dna_phred42 t32{'C'_rna4};
    [[maybe_unused]] constexpr qualified_gapped_dna_phred42 t33{'$'_phred42};
    [[maybe_unused]] constexpr qualified_gapped_dna_phred42 t34{seqan3::gap{}};
    [[maybe_unused]] constexpr qualified_gapped_dna_phred42 t35{seqan3::gapped<seqan3::dna4>('C'_dna4)};
    [[maybe_unused]] constexpr qualified_gapped_dna_phred42 t36{seqan3::gapped<seqan3::dna4>(seqan3::gap{})};
    [[maybe_unused]] constexpr qualified_gapped_dna_phred42 t37{seqan3::gap{}, '$'_phred42};

    [[maybe_unused]] constexpr gapped_qualified_dna_phred42 t41{'C'_dna4};
    [[maybe_unused]] constexpr gapped_qualified_dna_phred42 t42{'C'_rna4};
    [[maybe_unused]] constexpr gapped_qualified_dna_phred42 t43{'$'_phred42};
    [[maybe_unused]] constexpr gapped_qualified_dna_phred42 t44{seqan3::gap{}};
    [[maybe_unused]] constexpr gapped_qualified_dna_phred42 t45{qualified_dna_phred42{'C'_dna4, '!'_phred42}};

    [[maybe_unused]] constexpr qualified_qualified_gapped_dna_phred42_phred42 t51{'C'_dna4};
    [[maybe_unused]] constexpr qualified_qualified_gapped_dna_phred42_phred42 t52{'C'_rna4};
    [[maybe_unused]] constexpr qualified_qualified_gapped_dna_phred42_phred42 t53{'$'_phred42};
    [[maybe_unused]] constexpr qualified_qualified_gapped_dna_phred42_phred42 t54{seqan3::gap{}};
    [[maybe_unused]] constexpr qualified_qualified_gapped_dna_phred42_phred42 t55{
        seqan3::gapped<seqan3::dna4>('C'_dna4)};
    [[maybe_unused]] constexpr qualified_qualified_gapped_dna_phred42_phred42 t56{
        seqan3::gapped<seqan3::dna4>(seqan3::gap{})};

    [[maybe_unused]] constexpr gapped_alphabet_variant_dna_phred42 t61{'C'_dna4};
    [[maybe_unused]] constexpr gapped_alphabet_variant_dna_phred42 t62{'C'_rna4};
    [[maybe_unused]] constexpr gapped_alphabet_variant_dna_phred42 t63{'$'_phred42};
    [[maybe_unused]] constexpr gapped_alphabet_variant_dna_phred42 t64{seqan3::gap{}};
    [[maybe_unused]] constexpr gapped_alphabet_variant_dna_phred42 t65{qualified_dna_phred42{'C'_dna4, '!'_phred42}};
}

TEST(composite, custom_assignment)
{
    qualified_dna_phred42 t11{};
    qualified_dna_phred42 t12{'C'_dna4};
    qualified_dna_phred42 t13{'C'_dna4, '$'_phred42};
    t11 = 'C'_dna4;
    EXPECT_EQ(t11, t12);
    t11 = 'C'_rna4;
    EXPECT_EQ(t11, t12);
    t11 = '$'_phred42;
    EXPECT_EQ(t11, t13);
    // t11 = '$'_phred63; // does not work because of explicit conversion

    seqan3::qualified<seqan3::aa27, seqan3::phred63> t20{'K'_aa27, seqan3::phred63{}};
    seqan3::qualified<seqan3::aa27, seqan3::phred63> t21{};
    seqan3::qualified<seqan3::aa27, seqan3::phred63> t22{'K'_aa27, '$'_phred63};
    t21 = 'K'_aa27;
    EXPECT_EQ(t20, t21);
    t21 = '$'_phred63;
    EXPECT_EQ(t21, t22);

    qualified_gapped_dna_phred42 t31{};
    qualified_gapped_dna_phred42 t32{'C'_dna4};
    qualified_gapped_dna_phred42 t33{'C'_dna4, '$'_phred42};
    qualified_gapped_dna_phred42 t34{seqan3::gap{}, '$'_phred42};
    t31 = 'C'_dna4;
    EXPECT_EQ(t31, t32);
    t31 = 'C'_rna4;
    EXPECT_EQ(t31, t32);
    t31 = '$'_phred42;
    EXPECT_EQ(t31, t33);
    t31 = seqan3::gap{};
    EXPECT_EQ(t31, t34);
    t31 = seqan3::gapped<seqan3::dna4>('C'_dna4);
    EXPECT_EQ(t31, t33);
    t31 = seqan3::gapped<seqan3::dna4>(seqan3::gap{});
    EXPECT_EQ(t31, t34);

    gapped_qualified_dna_phred42 t41{};
    gapped_qualified_dna_phred42 t42{'C'_dna4};
    // gapped_qualified_dna_phred42 t43{qualified_dna_phred42{'C'_dna4, '$'_phred42}};
    gapped_qualified_dna_phred42 t44{seqan3::gap{}};
    gapped_qualified_dna_phred42 t45{qualified_dna_phred42{'C'_dna4, '!'_phred42}};
    t41 = 'C'_dna4;
    EXPECT_EQ(t41, t42);
    t41 = 'C'_rna4;
    EXPECT_EQ(t41, t42);
    t41 = '$'_phred42;
    // EXPECT_EQ(t41, t43); should work intuitively but does not because on assignment the seqan3::qualified object is defaulted
    t41 = seqan3::gap{};
    EXPECT_EQ(t41, t44);
    t41 = qualified_dna_phred42{'C'_dna4, '!'_phred42};
    EXPECT_EQ(t41, t45);

    qualified_qualified_gapped_dna_phred42_phred42 t51{};
    qualified_qualified_gapped_dna_phred42_phred42 t52{'C'_dna4};
    qualified_qualified_gapped_dna_phred42_phred42 t53{qualified_gapped_dna_phred42{'C'_dna4, '!'_phred42},
                                                       '$'_phred42};
    qualified_qualified_gapped_dna_phred42_phred42 t54{qualified_gapped_dna_phred42{seqan3::gap{}, '!'_phred42},
                                                       '$'_phred42};
    t51 = 'C'_dna4;
    EXPECT_EQ(t51, t52);
    t51 = 'C'_rna4;
    EXPECT_EQ(t51, t52);
    t51 = '$'_phred42;
    EXPECT_EQ(t51, t53);
    t51 = seqan3::gap{};
    EXPECT_EQ(t51, t54);
    t51 = seqan3::gapped<seqan3::dna4>('C'_dna4);
    EXPECT_EQ(t51, t53);
    t51 = seqan3::gapped<seqan3::dna4>(seqan3::gap{});
    EXPECT_EQ(t51, t54);

    gapped_alphabet_variant_dna_phred42 t61{};
    gapped_alphabet_variant_dna_phred42 t62{'C'_dna4};
    gapped_alphabet_variant_dna_phred42 t63{'$'_phred42};
    gapped_alphabet_variant_dna_phred42 t64{seqan3::gap{}};
    gapped_alphabet_variant_dna_phred42 t65{qualified_dna_phred42{'C'_dna4, '!'_phred42}};
    t61 = 'C'_dna4;
    EXPECT_EQ(t61, t62);
    t61 = 'C'_rna4;
    EXPECT_EQ(t61, t62);
    t61 = '$'_phred42;
    EXPECT_EQ(t61, t63);
    t61 = seqan3::gap{};
    EXPECT_EQ(t61, t64);
    t61 = qualified_dna_phred42{'C'_dna4, '!'_phred42};
    EXPECT_EQ(t61, t65);
}

constexpr bool do_assignment()
{
    qualified_dna_phred42 t11{};
    t11 = 'C'_dna4;
    t11 = 'C'_rna4;
    t11 = '$'_phred42;
    // t11 = '$'_phred63; // does not work because of explicit conversion

    seqan3::qualified<seqan3::aa27, seqan3::phred63> t21{};
    t21 = 'K'_aa27;
    t21 = '$'_phred63;

    qualified_gapped_dna_phred42 t31{};
    t31 = 'C'_dna4;
    t31 = 'C'_rna4;
    t31 = '$'_phred42;
    t31 = seqan3::gap{};
    t31 = seqan3::gapped<seqan3::dna4>('C'_dna4);
    t31 = seqan3::gapped<seqan3::dna4>(seqan3::gap{});

    gapped_qualified_dna_phred42 t41{};
    t41 = 'C'_dna4;
    t41 = 'C'_rna4;
    t41 = '$'_phred42;
    t41 = seqan3::gap{};
    t41 = qualified_dna_phred42{'C'_dna4, '!'_phred42};

    qualified_qualified_gapped_dna_phred42_phred42 t51{};
    t51 = 'C'_dna4;
    t51 = 'C'_rna4;
    t51 = '$'_phred42;
    t51 = seqan3::gap{};
    t51 = seqan3::gapped<seqan3::dna4>('C'_dna4);
    t51 = seqan3::gapped<seqan3::dna4>(seqan3::gap{});

    gapped_alphabet_variant_dna_phred42 t61{};
    t61 = 'C'_rna4;
    t61 = '$'_phred42;
    t61 = seqan3::gap{};
    t61 = qualified_dna_phred42{'C'_dna4, '!'_phred42};

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

    qualified_dna_phred42 t11{'C'_dna4, '$'_phred42};
    EXPECT_EQ(t11, 'C'_dna4);
    EXPECT_EQ(t11, 'C'_rna4);
    EXPECT_EQ(t11, '$'_phred42);
    EXPECT_LT(t11, 'G'_dna4);
    EXPECT_LT(t11, 'G'_rna4);
    EXPECT_LT(t11, '%'_phred42);

    EXPECT_EQ('C'_dna4, t11);
    EXPECT_EQ('C'_rna4, t11);
    EXPECT_EQ('$'_phred42, t11);
    EXPECT_LT('A'_dna4, t11);
    EXPECT_LT('A'_rna4, t11);
    EXPECT_LT('#'_phred42, t11);

    seqan3::qualified<seqan3::aa27, seqan3::phred63> t21{'K'_aa27, '$'_phred63};
    EXPECT_EQ(t21, 'K'_aa27);
    EXPECT_EQ(t21, '$'_phred63);
    EXPECT_LT(t21, 'L'_aa27);
    EXPECT_LT(t21, '%'_phred63);

    EXPECT_EQ('K'_aa27, t21);
    EXPECT_EQ('$'_phred63, t21);
    EXPECT_LT('C'_aa27, t21);
    EXPECT_LT('#'_phred63, t21);

    qualified_gapped_dna_phred42 t31{'C'_dna4, '$'_phred42};
    EXPECT_EQ(t31, 'C'_dna4);
    EXPECT_EQ(t31, 'C'_rna4);
    EXPECT_EQ(t31, '$'_phred42);
    EXPECT_NE(t31, seqan3::gap{});
    EXPECT_EQ(t31, seqan3::gapped<seqan3::dna4>('C'_dna4));
    EXPECT_LT(t31, 'G'_dna4); // *
    EXPECT_LT(t31, 'G'_rna4); // *
    EXPECT_LT(t31, '%'_phred42);
    EXPECT_LT(t31, seqan3::gap{}); // *
    EXPECT_LT(t31, seqan3::gapped<seqan3::dna4>('G'_dna4));

    EXPECT_EQ('C'_dna4, t31);
    EXPECT_EQ('C'_rna4, t31);
    EXPECT_EQ('$'_phred42, t31);
    EXPECT_NE(seqan3::gap{}, t31);
    EXPECT_EQ(seqan3::gapped<seqan3::dna4>('C'_dna4), t31);
    EXPECT_LT('A'_dna4, t31); // *
    EXPECT_LT('A'_rna4, t31); // *
    EXPECT_LT('#'_phred42, t31);
    EXPECT_GT(seqan3::gap{}, t31); // *
    EXPECT_LT(seqan3::gapped<seqan3::dna4>('A'_dna4), t31);

    gapped_qualified_dna_phred42 t41{qualified_dna_phred42{'C'_dna4, '$'_phred42}};
    EXPECT_EQ(t41, 'C'_dna4);
    EXPECT_EQ(t41, 'C'_rna4);
    EXPECT_EQ(t41, '$'_phred42);
    EXPECT_NE(t41, seqan3::gap{});
    EXPECT_EQ(t41, (qualified_dna_phred42{'C'_dna4, '$'_phred42}));
    EXPECT_EQ(t41, (gapped_qualified_dna_phred42{qualified_dna_phred42{'C'_dna4, '$'_phred42}}));
    //     EXPECT_LT(t41, 'G'_dna4); // not supposed to work
    //     EXPECT_LT(t41, 'G'_rna4); // not supposed to work
    //     EXPECT_LT(t41, '%'_phred42); // would never be LT, because seqan3::dna4 part of tuple defaulted to 'A' on RHS
    EXPECT_LT(t41, seqan3::gap{});                                  // *
    EXPECT_LT(t41, (qualified_dna_phred42{'G'_dna4, '#'_phred42})); // *
    EXPECT_LT(t41, (gapped_qualified_dna_phred42{qualified_dna_phred42{'G'_dna4, '#'_phred42}}));

    EXPECT_EQ('C'_dna4, t41);
    EXPECT_EQ('C'_rna4, t41);
    EXPECT_EQ('$'_phred42, t41);
    EXPECT_EQ((qualified_dna_phred42{'C'_dna4, '$'_phred42}), t41);
    EXPECT_NE(seqan3::gap{}, t41);
    //     EXPECT_LT('A'_dna4, t41); // not supposed to work
    //     EXPECT_LT('A'_rna4, t41); // not supposed to work
    //     EXPECT_LT('#'_phred42, t41); // not supposed to work
    EXPECT_LT((qualified_dna_phred42{'A'_dna4, '#'_phred42}), t41); // *
    EXPECT_GT(seqan3::gap{}, t41);                                  // *

    qualified_qualified_gapped_dna_phred42_phred42 t51{qualified_gapped_dna_phred42{'C'_dna4, '$'_phred42}};
    EXPECT_EQ(t51, 'C'_dna4);
    EXPECT_EQ(t51, 'C'_rna4);
    EXPECT_NE(t51, seqan3::gap{});
    EXPECT_EQ(t51, seqan3::gapped<seqan3::dna4>('C'_dna4));
    EXPECT_EQ(t51, '!'_phred42); // "outer" Phred score element
    EXPECT_EQ(t51, (qualified_gapped_dna_phred42{'C'_dna4, '$'_phred42}));
    //     EXPECT_LT(t51, 'G'_dna4); // not supposed to work
    //     EXPECT_LT(t51, 'G'_rna4); // not supposed to work
    //     EXPECT_LT(t51, seqan3::gap{}); // not supposed to work
    //     EXPECT_LT(t51, seqan3::gapped<seqan3::dna4>('G'_dna4)); // not supposed to work
    EXPECT_LT(t51, '"'_phred42);
    EXPECT_LT(t51, (qualified_gapped_dna_phred42{'C'_dna4, '%'_phred42}));

    EXPECT_EQ('C'_dna4, t51);
    EXPECT_EQ('C'_rna4, t51);
    EXPECT_NE(seqan3::gap{}, t51);
    EXPECT_EQ(seqan3::gapped<seqan3::dna4>('C'_dna4), t51);
    EXPECT_EQ('!'_phred42, t51);
    EXPECT_EQ((qualified_gapped_dna_phred42{'C'_dna4, '$'_phred42}), t51);
    //     EXPECT_LT('A'_dna4, t51); // not supposed to work
    //     EXPECT_LT('A'_rna4, t51); // not supposed to work
    //     EXPECT_GT(seqan3::gap{}, t51); // not supposed to work
    //     EXPECT_LT(seqan3::gapped<seqan3::dna4>('A'_dna4), t51); // not supposed to work
    EXPECT_GT('"'_phred42, t51);
    EXPECT_GT((qualified_gapped_dna_phred42{'C'_dna4, '%'_phred42}), t51);

    gapped_alphabet_variant_dna_phred42 t61{'C'_rna4};
    EXPECT_EQ(t61, 'C'_rna4);
    EXPECT_EQ(t61, 'C'_dna4);
    EXPECT_NE(t61, seqan3::gap{});
    EXPECT_NE(t61, '!'_phred42);
    EXPECT_LT(t61, 'G'_rna4);      // *
    EXPECT_LT(t61, 'G'_dna4);      // *
    EXPECT_LT(t61, seqan3::gap{}); // *
    EXPECT_LT(t61, '"'_phred42);   // *

    EXPECT_EQ('C'_rna4, t61);
    EXPECT_EQ('C'_dna4, t61);
    EXPECT_NE(seqan3::gap{}, t61);
    EXPECT_NE('!'_phred42, t61);
    EXPECT_LT('A'_rna4, t61);      // *
    EXPECT_LT('A'_dna4, t61);      // *
    EXPECT_GT(seqan3::gap{}, t61); // *
    EXPECT_GT('!'_phred42, t61);   // *
}

TEST(composite, get)
{
    using seqan3::get;

    qualified_qualified_gapped_dna_phred42_phred42 t51{qualified_gapped_dna_phred42{'C'_dna4, '$'_phred42}};
    EXPECT_EQ(get<0>(t51), 'C'_dna4);
    EXPECT_EQ(get<0>(get<0>(t51)), 'C'_dna4);

    EXPECT_EQ(get<0>(t51), 'C'_rna4);
    EXPECT_EQ(get<0>(get<0>(t51)), 'C'_rna4);

    EXPECT_NE(get<0>(t51), seqan3::gap{});
    EXPECT_NE(get<0>(get<0>(t51)), seqan3::gap{});

    EXPECT_EQ(get<0>(t51), seqan3::gapped<seqan3::dna4>('C'_dna4));
    EXPECT_EQ(get<0>(get<0>(t51)), seqan3::gapped<seqan3::dna4>('C'_dna4));

    EXPECT_NE(get<0>(t51), '!'_phred42);
}
