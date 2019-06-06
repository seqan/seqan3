// -----------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2019, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2019, Knut Reinert & MPI für molekulare Genetik
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at: https://github.com/seqan/seqan3/blob/master/LICENSE.md
// -----------------------------------------------------------------------------------------------------

#include <type_traits>

#include <gtest/gtest.h>

#include <seqan3/core/type_traits/concept.hpp>
#include <seqan3/core/type_traits/template_inspection.hpp>
#include <seqan3/core/type_list.hpp>

using namespace seqan3;

TEST(template_inspect, concept_check)
{
    using tl = type_list<int, char, double>;

    EXPECT_FALSE((TransformationTrait<detail::transfer_template_args_onto<int, std::tuple>>));
    EXPECT_TRUE((TransformationTrait<detail::transfer_template_args_onto<tl, std::tuple>>));

    EXPECT_TRUE((UnaryTypeTrait<detail::is_type_specialisation_of<int, type_list>>));
}

TEST(template_inspect, transfer_template_args_onto_t)
{
    using tl = type_list<int, char, double>;
    using t = detail::transfer_template_args_onto<tl, std::tuple>::type;
    EXPECT_TRUE((std::is_same_v<t, std::tuple<int, char, double>>));

    // shortcut
    using t = detail::transfer_template_args_onto_t<tl, std::tuple>;
    EXPECT_TRUE((std::is_same_v<t, std::tuple<int, char, double>>));
}

TEST(template_inspect, is_type_specialisation_of)
{
    using tl = type_list<int, char, double>;
    EXPECT_TRUE((detail::is_type_specialisation_of<tl, type_list>::value));
    EXPECT_FALSE((detail::is_type_specialisation_of<int, type_list>::value));
}

TEST(template_inspect, is_type_specialisation_of_v)
{
    using tl = type_list<int, char, double>;
    EXPECT_TRUE((detail::is_type_specialisation_of_v<tl, type_list>));
    EXPECT_FALSE((detail::is_type_specialisation_of_v<int, type_list>));
}

template <std::Integral t>
struct constraint_bar
{};

TEST(template_inspect, is_type_specialisation_of_with_ill_formed_type)
{
    EXPECT_FALSE((detail::is_type_specialisation_of<std::tuple<float>, constraint_bar>::value));
}

template <int i, char c>
struct t1 {};

template <int _i, char _c>
struct t2
{
    static constexpr auto i = _i;
    static constexpr auto c = _c;
};

enum struct e1
{
    foo
};

template <e1 v>
struct foo
{};

enum struct e2
{
    bar
};

template <e2 v>
struct bar
{};

template <e2 v>
struct bar2
{};

TEST(template_inspect, transfer_template_vargs_onto_enum)
{
    EXPECT_TRUE((std::is_same_v<detail::transformation_trait_or_t<
                                    detail::transfer_template_vargs_onto<bar<e2::bar>, foo>,
                                    void>,
                                void>));

    using ta2 = detail::transfer_template_vargs_onto<bar<e2::bar>, bar>::type;
    EXPECT_TRUE((std::is_same_v<ta2, bar<e2::bar>>));

    using ta3 = detail::transfer_template_vargs_onto<bar<e2::bar>, bar2>::type;
    EXPECT_TRUE((std::is_same_v<ta3, bar2<e2::bar>>));
}

TEST(template_inspect, transfer_template_vargs_onto_t)
{
    using tl = t1<1, 'a'>;
    using ta = detail::transfer_template_vargs_onto<tl, t2>::type;
    EXPECT_EQ(1,   ta::i);
    EXPECT_EQ('a', ta::c);

    // shortcut
    using ta2 = detail::transfer_template_vargs_onto_t<tl, t2>;
    EXPECT_EQ(1,   ta2::i);
    EXPECT_EQ('a', ta2::c);
}

TEST(template_inspect, is_value_specialisation_of)
{
    using tl = t1<1, 'a'>;

    EXPECT_TRUE((detail::is_value_specialisation_of<tl, t1>::value));
    EXPECT_FALSE((detail::is_value_specialisation_of<int, t1>::value));
}

TEST(template_inspect, is_value_specialisation_of_v)
{
    using tl = t1<1, 'a'>;

    EXPECT_TRUE((detail::is_value_specialisation_of_v<tl, t1>));
    EXPECT_FALSE((detail::is_value_specialisation_of_v<int, t1>));
}

template <int varg>
    requires  0 <= varg && varg <= 2
struct constraint_vbar
{};

template <auto ...vargs>
struct vargs_foo
{};

TEST(template_inspect, is_type_specialisation_of_with_ill_formed_non_type_template)
{
    EXPECT_FALSE((detail::is_value_specialisation_of_v<vargs_foo<5>, constraint_vbar>));
}
