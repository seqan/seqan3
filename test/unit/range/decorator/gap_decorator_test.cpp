// -----------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2020, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2020, Knut Reinert & MPI für molekulare Genetik
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at: https://github.com/seqan/seqan3/blob/master/LICENSE.md
// -----------------------------------------------------------------------------------------------------

#include <seqan3/std/ranges>
#include <vector>

#include <gtest/gtest.h>

#include <range/v3/view/filter.hpp>

#include <seqan3/alignment/aligned_sequence/aligned_sequence_concept.hpp>
#include <seqan3/alphabet/nucleotide/dna4.hpp>
#include <seqan3/range/decorator/gap_decorator.hpp>
#include <seqan3/range/views/enforce_random_access.hpp>
#include <seqan3/range/views/to_char.hpp>

#include "../iterator_test_template.hpp"
#include "../../alignment/aligned_sequence_test_template.hpp"

using decorator_t = seqan3::gap_decorator<std::vector<seqan3::dna4> const &>;

std::vector<seqan3::dna4> const dummy_obj{}; // dummy lvalue for type declaration of views
using decorator_t2 = seqan3::gap_decorator<
                         decltype(std::ranges::subrange<decltype(dummy_obj.begin()),
                                                        decltype(dummy_obj.begin())>{dummy_obj.begin(),
                                                                                     dummy_obj.end()})>;


using test_types = ::testing::Types<decorator_t, decorator_t2>;

// ---------------------------------------------------------------------------------------------------------------------
// test templates
// ---------------------------------------------------------------------------------------------------------------------

template <typename inner_type_>
class aligned_sequence_<seqan3::gap_decorator<inner_type_>> : public ::testing::Test
{
public:
    // Initialiser function is needed for the typed test because the gapped_decorator
    // will be initialised differently than the naive vector<seqan3::gapped<seqan3::dna>>.
    void initialise_typed_test_container(decorator_t & container_, seqan3::dna4_vector const & target)
    {
        container_ = target;
    }

    // Initialiser function is needed for the typed test because the gapped_decorator
    // will be initialised differently than the naive vector<seqan3::gapped<seqan3::dna>>.
    void initialise_typed_test_container(decorator_t2 & container_, seqan3::dna4_vector const & target)
    {
        container_ = std::ranges::subrange<decltype(target.begin()),
                                           decltype(target.end())>{target.begin(), target.end()};
    }
};

INSTANTIATE_TYPED_TEST_SUITE_P(gap_decorator, aligned_sequence_, test_types, );

template <>
struct iterator_fixture<std::ranges::iterator_t<decorator_t>> : public ::testing::Test
{
    using iterator_tag = std::bidirectional_iterator_tag;
    static constexpr bool const_iterable = true;

    std::vector<seqan3::dna4> const vec{"ACTGACTG"_dna4};

    template <typename t>
    static void initialise_with_gaps(t & v)
    {
        insert_gap(v, std::next(v.begin(), 5), 4);
        insert_gap(v, std::next(v.begin(), 2));
        insert_gap(v, v.end(), 3);
        insert_gap(v, v.begin(), 5);
    }

    std::vector<seqan3::gapped<seqan3::dna4>> expected_range = [this] ()
    {
        std::vector<seqan3::gapped<seqan3::dna4>> tmp{};
        assign_unaligned(tmp, vec);
        initialise_with_gaps(tmp);
        return tmp;
    }();

    decorator_t test_range = [this] ()
    {
        decorator_t tmp{vec};
        initialise_with_gaps(tmp);
        return tmp;
    }();
};

struct random_access_iterator_test
{};

template <>
struct iterator_fixture<random_access_iterator_test> : iterator_fixture<std::ranges::iterator_t<decorator_t>>
{
    using base_t = iterator_fixture<std::ranges::iterator_t<decorator_t>>;

    using iterator_tag = std::random_access_iterator_tag;
    static constexpr bool const_iterable = true;

    decltype(base_t::test_range
             | seqan3::views::enforce_random_access) test_range = base_t::test_range
                                                                | seqan3::views::enforce_random_access;
};

INSTANTIATE_TYPED_TEST_SUITE_P(gap_decorator_iterator, iterator_fixture, std::ranges::iterator_t<decorator_t>, );
INSTANTIATE_TYPED_TEST_SUITE_P(gap_decorator_iterator_random_access, iterator_fixture, random_access_iterator_test, );

// ---------------------------------------------------------------------------------------------------------------------
// typed test
// ---------------------------------------------------------------------------------------------------------------------

template <typename t>
class gap_decorator_f : public ::testing::Test {};

TYPED_TEST_SUITE(gap_decorator_f, test_types, );

// concept checks
TYPED_TEST(gap_decorator_f, concept_checks)
{
    EXPECT_TRUE((std::ranges::bidirectional_range<TypeParam>));
    EXPECT_TRUE((std::ranges::bidirectional_range<TypeParam const>));

    EXPECT_FALSE((std::ranges::enable_view<TypeParam>));
    EXPECT_FALSE((std::ranges::enable_view<TypeParam &>));
    EXPECT_FALSE((ranges::enable_view<TypeParam>));
    EXPECT_FALSE((ranges::enable_view<TypeParam &>));

    EXPECT_FALSE((std::ranges::view<TypeParam>));

    EXPECT_TRUE((seqan3::aligned_sequence<TypeParam>));
}

TYPED_TEST(gap_decorator_f, construction_general)
{
    // default
    {
        [[maybe_unused]] TypeParam dec{};
        [[maybe_unused]] TypeParam dec3{};
    }

    // copy
    {
        [[maybe_unused]] TypeParam dec{};
        [[maybe_unused]] TypeParam dec2(dec);

        [[maybe_unused]] TypeParam const dec_const{};
        [[maybe_unused]] TypeParam const dec2_const(dec_const);
    }

    // move
    {
        [[maybe_unused]] TypeParam dec{};
        [[maybe_unused]] TypeParam dec2(std::move(dec));

        // const
        [[maybe_unused]] TypeParam const dec_const{};
        [[maybe_unused]] TypeParam const dec2_const(std::move(dec_const));
    }

    // copy assignment
    {
        [[maybe_unused]] TypeParam dec{};
        [[maybe_unused]] TypeParam dec2;
        dec2 = dec;
    }

    // move assignment
    {
        [[maybe_unused]] TypeParam dec{};
        [[maybe_unused]] TypeParam dec2;
        dec2 = std::move(dec);
    }
}

// ---------------------------------------------------------------------------------------------------------------------
// general test with automatic type deduction
// ---------------------------------------------------------------------------------------------------------------------

TEST(gap_decorator, construction_from_ungapped_sequence)
{
    {
        std::vector<seqan3::dna4> v{"ACTG"_dna4};
        std::vector<seqan3::dna4> const v_const{"ACTG"_dna4};

        // non-const version
        seqan3::gap_decorator dec{v};
        EXPECT_EQ('A'_dna4, dec[0]);
        EXPECT_EQ('C'_dna4, dec[1]);

        // const version
        seqan3::gap_decorator const dec2(v_const);
        EXPECT_EQ('A'_dna4, dec2[0]);
        EXPECT_EQ('C'_dna4, dec2[1]);
    }

    {
        std::vector<seqan3::dna4> v{"ACTG"_dna4};
        std::vector<seqan3::dna4> const v_const{"ACTG"_dna4};

        // non-const version
        seqan3::gap_decorator dec = v;
        EXPECT_EQ('A'_dna4, dec[0]);
        EXPECT_EQ('C'_dna4, dec[1]);

        // const version
        seqan3::gap_decorator const dec2 = v_const;
        EXPECT_EQ('A'_dna4, dec2[0]);
        EXPECT_EQ('C'_dna4, dec2[1]);
    }
}

TEST(gap_decorator, assignment_from_ungapped_sequence)
{
    {
        std::vector<seqan3::dna4> v{"TT"_dna4};
        std::vector<seqan3::dna4> v2{"ACTG"_dna4};
        std::vector<seqan3::dna4> const v_const{"TGCC"_dna4};

        seqan3::gap_decorator dec{v};
        dec = v2;
        EXPECT_EQ('A'_dna4, dec[0]);
        EXPECT_EQ('C'_dna4, dec[1]);

        // const vec
        dec = v_const;
        EXPECT_EQ('T'_dna4, dec[0]);
        EXPECT_EQ('G'_dna4, dec[1]);

        // reassignment after adding gaps
        EXPECT_EQ(dec.size(), v2.size());
        insert_gap(dec, begin(dec), 2);
        EXPECT_EQ(dec.size(), v2.size() + 2);
        dec = v2;
        EXPECT_EQ(dec.size(), v2.size());
        EXPECT_EQ('A'_dna4, dec[0]);
        EXPECT_EQ('C'_dna4, dec[1]);
    }
}

TEST(gap_decorator, comparison)
{
    std::vector<seqan3::dna4> v{"ACTG"_dna4};

    seqan3::gap_decorator dec{v};
    seqan3::gap_decorator dec2{v};

    EXPECT_EQ(dec, dec2);
    EXPECT_LE(dec, dec2);
    EXPECT_GE(dec, dec2);

    insert_gap(dec, dec.end(), 2);

    EXPECT_NE(dec, dec2);
    EXPECT_LT(dec2, dec); // dec2 is prefix of dec
    EXPECT_LE(dec2, dec); // dec2 is prefix of dec
    EXPECT_GT(dec, dec2); // dec2 is prefix of dec
    EXPECT_GE(dec, dec2); // dec2 is prefix of dec

    insert_gap(dec2, dec2.end(), 2);
    insert_gap(dec2, dec2.begin(), 1);

    EXPECT_NE(dec, dec2); // ACTG-- vs -ACTG--
    EXPECT_GT(dec2, dec);
    EXPECT_GE(dec2, dec);
    EXPECT_LT(dec, dec2);
    EXPECT_LE(dec, dec2);

    std::vector<seqan3::dna4> v2{"TCTG"_dna4};
    seqan3::gap_decorator decNE{v2};
    EXPECT_NE(dec, decNE);
}

TEST(gap_decorator, begin_and_end)
{
    std::vector<seqan3::dna4> v{"ACTG"_dna4};

    seqan3::gap_decorator dec{v};
    seqan3::gap_decorator const dec_const{v};

    EXPECT_EQ(*dec.begin(), 'A'_dna4);
    EXPECT_EQ(*dec.cbegin(), 'A'_dna4);
    EXPECT_EQ(*dec_const.begin(), 'A'_dna4);
    EXPECT_EQ(*dec_const.cbegin(), 'A'_dna4);

    [[maybe_unused]] auto end = dec.end();
    [[maybe_unused]] auto end_const = dec.cend();
}

TEST(gap_decorator, decorator_on_views)
{
    std::vector<seqan3::dna4> v{"ACTG"_dna4};

    auto sub = std::ranges::subrange<decltype(v.begin()), decltype(v.begin())>{v.begin()+1, v.begin()+3};
    seqan3::gap_decorator dec{sub};

    EXPECT_EQ(dec.size(), 2u);
    EXPECT_EQ(*dec.begin(), 'C'_dna4);
    EXPECT_EQ(dec[1], 'T'_dna4);

    auto it = insert_gap(dec, std::next(dec.begin(), 1), 2);

    EXPECT_EQ(dec.size(), 4u);
    EXPECT_EQ(*dec.begin(), 'C'_dna4);
    EXPECT_EQ(*(std::next(dec.begin(), 1)), seqan3::gap{});
    EXPECT_EQ(*it, seqan3::gap{});

    // auto v_char = v | seqan3::views::to_char;
    seqan3::gap_decorator dec2{v | seqan3::views::to_char};
    EXPECT_EQ(dec2.size(), 4u);
    EXPECT_EQ(*dec2.begin(), 'A');
    EXPECT_EQ(*++dec2.begin(), 'C');

    auto dec3 = dec | std::views::filter([] (auto chr) { return chr != seqan3::gap{}; });
    EXPECT_EQ(*dec3.begin(), 'C'_dna4);
    EXPECT_EQ(*(std::next(dec3.begin())), 'T'_dna4);
}
