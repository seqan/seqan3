// SPDX-FileCopyrightText: 2006-2025 Knut Reinert & Freie Universität Berlin
// SPDX-FileCopyrightText: 2016-2025 Knut Reinert & MPI für molekulare Genetik
// SPDX-License-Identifier: BSD-3-Clause

#include <gtest/gtest.h>

#include <type_traits>

#include <seqan3/core/detail/template_inspection.hpp>
#include <seqan3/utility/type_traits/concept.hpp>

template <typename... args_t>
struct my_list;

template <std::integral t>
struct constraint_bar
{
    static_assert(std::same_as<t, int>);
};

TEST(template_inspect, concept_check)
{
    using tl = my_list<int, char, double>;

    EXPECT_FALSE((seqan3::transformation_trait<seqan3::detail::transfer_template_args_onto<int, std::tuple>>));
    EXPECT_TRUE((seqan3::transformation_trait<seqan3::detail::transfer_template_args_onto<tl, std::tuple>>));

    EXPECT_TRUE((seqan3::unary_type_trait<seqan3::detail::is_type_specialisation_of<int, my_list>>));
}

TEST(template_inspect, transfer_template_args_onto_t)
{
    using tl = my_list<int, char, double>;
    using t = seqan3::detail::transfer_template_args_onto<tl, std::tuple>::type;
    EXPECT_TRUE((std::is_same_v<t, std::tuple<int, char, double>>));

    // shortcut
    using t = seqan3::detail::transfer_template_args_onto_t<tl, std::tuple>;
    EXPECT_TRUE((std::is_same_v<t, std::tuple<int, char, double>>));
}

TEST(template_inspect, transfer_template_args_onto_with_constraint)
{
    // This test ensures that transfer_template_args_onto uses internally only
    // type declarations and not type instantiations

    using bar_char = constraint_bar<char>; // this is fine, even though it contains a static_assert
    using bar_int = constraint_bar<int>;

    // float does not fulfil integral constraint
    using bar_float_identity = seqan3::detail::transfer_template_args_onto<my_list<float>, constraint_bar>;
    EXPECT_FALSE(seqan3::transformation_trait<bar_float_identity>);

    // int fulfils integral constraint and static_assert
    using bar_int_identity = seqan3::detail::transfer_template_args_onto<my_list<int>, constraint_bar>;
    EXPECT_TRUE(seqan3::transformation_trait<bar_int_identity>);
    EXPECT_TRUE((std::is_same_v<typename bar_int_identity::type, bar_int>));

    // char fulfils integral constraint, but not static_assert
    using bar_char_identity = seqan3::detail::transfer_template_args_onto<my_list<char>, constraint_bar>;
    EXPECT_TRUE(seqan3::transformation_trait<bar_char_identity>);
    EXPECT_TRUE((std::is_same_v<typename bar_char_identity::type, bar_char>));
}

TEST(template_inspect, is_type_specialisation_of)
{
    using tl = my_list<int, char, double>;
    EXPECT_TRUE((seqan3::detail::is_type_specialisation_of<tl, my_list>::value));
    EXPECT_FALSE((seqan3::detail::is_type_specialisation_of<int, my_list>::value));
}

TEST(template_inspect, is_type_specialisation_of_v)
{
    using tl = my_list<int, char, double>;
    EXPECT_TRUE((seqan3::detail::is_type_specialisation_of_v<tl, my_list>));
    EXPECT_FALSE((seqan3::detail::is_type_specialisation_of_v<int, my_list>));
}

TEST(template_inspect, is_type_specialisation_with_constraint)
{
    // This test ensures that is_type_specialisation_of uses internally only
    // type declarations and not type instantiations

    using bar_char = constraint_bar<char>; // this is fine, even though it contains a static_assert
    using bar_int = constraint_bar<char>;

    EXPECT_FALSE((seqan3::detail::is_type_specialisation_of<std::tuple<float>, constraint_bar>::value));
    EXPECT_FALSE((seqan3::detail::is_type_specialisation_of<std::tuple<int>, constraint_bar>::value));

    EXPECT_TRUE((seqan3::detail::is_type_specialisation_of<bar_char, constraint_bar>::value));
    EXPECT_TRUE((seqan3::detail::is_type_specialisation_of<bar_int, constraint_bar>::value));
}

template <int i, char c>
struct t1
{};

template <int _i, char _c>
struct t2
{
    static constexpr auto i = _i;
    static constexpr auto c = _c;

    static_assert(c == 'a');
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
    bar,
    baz
};

template <e2 v>
struct bar
{};

template <e2 v>
struct bar2
{
    static_assert(v == e2::baz);
};

TEST(template_inspect, transfer_template_vargs_onto_enum)
{
    using foo_e2_bar = seqan3::detail::transfer_template_vargs_onto<bar<e2::bar>, foo>;
    EXPECT_TRUE((std::is_same_v<seqan3::detail::transformation_trait_or_t<foo_e2_bar, void>, void>));

    using ta2 = seqan3::detail::transfer_template_vargs_onto<bar<e2::bar>, bar>::type;
    EXPECT_TRUE((std::is_same_v<ta2, bar<e2::bar>>));

    // ensures that transfer_template_vargs_onto uses internally only type declarations and not type instantiations
    using ta3 = seqan3::detail::transfer_template_vargs_onto<bar<e2::bar>, bar2>::type;
    EXPECT_TRUE((std::is_same_v<ta3, bar2<e2::bar>>));
}

TEST(template_inspect, transfer_template_vargs_onto_t)
{
    using ta = seqan3::detail::transfer_template_vargs_onto<t1<1, 'a'>, t2>::type;
    EXPECT_EQ(1, ta::i);
    EXPECT_EQ('a', ta::c);

    // ensures that transfer_template_vargs_onto uses internally only type declarations and not type instantiations
    // 'b' does not fulfil the static_assert != 'a'
    using t2_identity = seqan3::detail::transfer_template_vargs_onto<t1<10, 'b'>, t2>;
    EXPECT_TRUE(seqan3::transformation_trait<t2_identity>);
    EXPECT_TRUE((std::is_same_v<typename t2_identity::type, t2<10, 'b'>>));
    // typename t2_identity::type{}; // instantiation of the type will static_assert

    // shortcut
    using ta2 = seqan3::detail::transfer_template_vargs_onto_t<t1<2, 'a'>, t2>;
    EXPECT_EQ(2, ta2::i);
    EXPECT_EQ('a', ta2::c);
}

TEST(template_inspect, is_value_specialisation_of)
{
    EXPECT_TRUE((seqan3::detail::is_value_specialisation_of<t1<1, 'a'>, t1>::value));
    EXPECT_FALSE((seqan3::detail::is_value_specialisation_of<int, t1>::value));
}

TEST(template_inspect, is_value_specialisation_of_v)
{
    EXPECT_TRUE((seqan3::detail::is_value_specialisation_of_v<t2<1, 'a'>, t2>));

    // ensures that is_value_specialisation_of_v uses internally only type declarations and not type instantiations
    EXPECT_TRUE((seqan3::detail::is_value_specialisation_of_v<t2<1, 'b'>, t2>));

    EXPECT_FALSE((seqan3::detail::is_value_specialisation_of_v<int, t1>));
}

template <int varg>
    requires (0 <= varg && varg <= 2)
struct constraint_vbar
{};

template <auto... vargs>
struct vargs_foo
{};

TEST(template_inspect, is_type_specialisation_of_with_ill_formed_non_type_template)
{
    EXPECT_FALSE((seqan3::detail::is_value_specialisation_of_v<vargs_foo<5>, constraint_vbar>));
}

TEST(template_inspect, template_specialisation_of)
{
    EXPECT_TRUE((seqan3::detail::template_specialisation_of<my_list<float>, my_list>));
    EXPECT_FALSE((seqan3::detail::template_specialisation_of<my_list<int>, std::tuple>));
}
