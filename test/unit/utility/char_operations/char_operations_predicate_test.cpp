// SPDX-FileCopyrightText: 2006-2025 Knut Reinert & Freie Universität Berlin
// SPDX-FileCopyrightText: 2016-2025 Knut Reinert & MPI für molekulare Genetik
// SPDX-License-Identifier: BSD-3-Clause

#include <gtest/gtest.h>

#include <string>

#include <seqan3/alphabet/aminoacid/aa27.hpp>
#include <seqan3/alphabet/concept.hpp>
#include <seqan3/alphabet/gap/gapped.hpp>
#include <seqan3/alphabet/nucleotide/dna4.hpp>
#include <seqan3/alphabet/nucleotide/rna5.hpp>
#include <seqan3/utility/char_operations/predicate.hpp>

using namespace std::literals;

template <char char_v>
struct foo : seqan3::detail::char_predicate_base<foo<char_v>>
{
    static inline std::string const msg{std::string{"foo_"} + std::string{char_v}};

    using base_t = seqan3::detail::char_predicate_base<foo<char_v>>;

    static constexpr auto data = []()
    {
        typename base_t::data_t d{};
        d[char_v] = true;
        return d;
    }();
};

template <char char_v>
inline constexpr foo<char_v> foo_v{};

struct bar
{
    template <typename alphabet_t>
    bool operator()([[maybe_unused]] alphabet_t c) const
    {
        return true;
    }
};

TEST(char_predicate, basic)
{
    foo<'a'> p{};
    EXPECT_TRUE(p('a'));
    EXPECT_FALSE(p('f'));
}

TEST(char_predicate, char_predicate_msg)
{
    EXPECT_EQ(foo<'o'>::msg, "foo_o"s);
}

TEST(char_predicate, concept)
{
    using seqan3::operator""_aa27;

    EXPECT_TRUE(seqan3::detail::char_predicate<decltype(seqan3::is_char<seqan3::to_char('A'_aa27)>)>);
    EXPECT_TRUE((seqan3::detail::char_predicate<decltype(seqan3::is_in_interval<'a', 'z'>)>));
    EXPECT_TRUE(seqan3::detail::char_predicate<decltype(seqan3::is_space)>);
    EXPECT_TRUE(seqan3::detail::char_predicate<decltype(seqan3::is_blank)>);
    EXPECT_TRUE(seqan3::detail::char_predicate<decltype(seqan3::is_graph)>);
    EXPECT_TRUE(seqan3::detail::char_predicate<decltype(seqan3::is_alpha)>);
    EXPECT_TRUE(seqan3::detail::char_predicate<decltype(seqan3::is_digit)>);
    EXPECT_TRUE(seqan3::detail::char_predicate<decltype(seqan3::is_alnum)>);
    // TODO(rrahn): Direct call in concept caused ICE.
    auto val = ((!seqan3::is_space || seqan3::is_alpha) || seqan3::is_digit);
    EXPECT_TRUE(seqan3::detail::char_predicate<decltype(val)>);
    EXPECT_TRUE(seqan3::detail::char_predicate<foo<' '>>);
    EXPECT_FALSE(seqan3::detail::char_predicate<bar>);
    EXPECT_FALSE(seqan3::detail::char_predicate<int>);
}

TEST(char_predicate, char_predicate_disjunction)
{
    using cond_t = seqan3::detail::char_predicate_disjunction<foo<'a'>, foo<'A'>, foo<'0'>>;
    EXPECT_TRUE(cond_t{}('a'));
    EXPECT_TRUE(cond_t{}('A'));
    EXPECT_TRUE(cond_t{}('0'));
    EXPECT_FALSE(cond_t{}('z'));
    EXPECT_FALSE(cond_t{}('!'));
    EXPECT_FALSE(cond_t{}('1'));

    constexpr auto p = foo_v<'a'> || foo_v<'A'> || foo_v<'0'>;
    EXPECT_TRUE(p('a'));
    EXPECT_TRUE(p('A'));
    EXPECT_TRUE(p('0'));
    EXPECT_FALSE(p('z'));
    EXPECT_FALSE(p('!'));
    EXPECT_FALSE(p('1'));
}

TEST(char_predicate, char_predicate_disjunction_msg)
{
    using or_t = seqan3::detail::char_predicate_disjunction<foo<'a'>, foo<'A'>, foo<'0'>>;
    EXPECT_EQ(or_t::msg, "(foo_a || foo_A || foo_0)"s);
}

TEST(char_predicate, is_not)
{
    using cond_t = seqan3::detail::char_predicate_negator<foo<'a'>>;
    EXPECT_FALSE(cond_t{}('a'));
    EXPECT_TRUE(cond_t{}('A'));
    EXPECT_TRUE(cond_t{}('0'));

    constexpr auto p = !foo_v<'a'>;
    EXPECT_FALSE(p('a'));
    EXPECT_TRUE(p('A'));
    EXPECT_TRUE(p('0'));
}

TEST(char_predicate, is_not_msg)
{
    using fn = decltype(!seqan3::is_alpha);
    EXPECT_EQ(fn::msg, "!(is_in_interval<'A', 'Z'> || is_in_interval<'a', 'z'>)"s);
}

TEST(char_predicate, is_in_interval)
{
    constexpr auto cond = seqan3::is_in_interval<'a', 'z'>;
    EXPECT_TRUE(cond('a'));
    EXPECT_TRUE(cond('k'));
    EXPECT_TRUE(cond('z'));
    EXPECT_FALSE(cond('A'));
    EXPECT_FALSE(cond('0'));
    EXPECT_FALSE(cond('!'));
}

TEST(char_predicate, is_in_interval_msg)
{
    EXPECT_EQ((seqan3::detail::is_in_interval_type<'a', 'z'>::msg), "is_in_interval<'a', 'z'>"s);
}

TEST(char_predicate, is_char)
{
    using seqan3::operator""_aa27;
    {
        constexpr auto cond = seqan3::is_char<'A'>;
        EXPECT_TRUE(cond('A'));
        EXPECT_FALSE(cond('x'));
    }

    {
        constexpr auto cond = seqan3::is_char<seqan3::to_char('A'_aa27)>;
        EXPECT_TRUE(cond('A'));
        EXPECT_FALSE(cond('z'));
    }
}

TEST(char_predicate, is_char_msg)
{
    using seqan3::operator""_dna4;
    EXPECT_EQ((seqan3::is_char<seqan3::to_char('A'_dna4)>.msg), "is_char<'A'>"s);
    EXPECT_EQ((seqan3::is_char<'\t'>.msg), "is_char<'\t'>"s);
}

TEST(char_predicate, is_cntrl)
{
    EXPECT_TRUE(seqan3::is_cntrl('\0'));
    EXPECT_TRUE(seqan3::is_cntrl(static_cast<char>(31)));
    EXPECT_TRUE(seqan3::is_cntrl(static_cast<char>(127)));
    EXPECT_TRUE(seqan3::is_cntrl('\t'));
    EXPECT_FALSE(seqan3::is_cntrl('A'));
}

TEST(char_predicate, is_print)
{
    EXPECT_FALSE(seqan3::is_print('\0'));
    EXPECT_FALSE(seqan3::is_print(static_cast<char>(31)));
    EXPECT_FALSE(seqan3::is_print(static_cast<char>(127)));
    EXPECT_TRUE(seqan3::is_print(' '));
    EXPECT_TRUE(seqan3::is_print('A'));
    EXPECT_TRUE(seqan3::is_print('~'));
}

TEST(char_predicate, is_print_msg)
{
    EXPECT_EQ((seqan3::is_print.message()), "is_in_interval<' ', '~'>"s);
}

TEST(char_predicate, is_blank)
{
    EXPECT_TRUE(seqan3::is_blank(' '));
    EXPECT_TRUE(seqan3::is_blank('\t'));
    EXPECT_FALSE(seqan3::is_blank('A'));
    EXPECT_FALSE(seqan3::is_blank('\n'));
}

TEST(char_predicate, is_blank_msg)
{
    EXPECT_EQ((seqan3::is_blank.message()), "(is_char<'\t'> || is_char<' '>)"s);
}

TEST(char_predicate, is_space)
{
    EXPECT_TRUE(seqan3::is_space('\n'));
    EXPECT_TRUE(seqan3::is_space('\r'));
    EXPECT_TRUE(seqan3::is_space('\f'));
    EXPECT_TRUE(seqan3::is_space('\v'));
    EXPECT_TRUE(seqan3::is_space('\t'));
    EXPECT_TRUE(seqan3::is_space(' '));
    EXPECT_FALSE(seqan3::is_space('0'));
    EXPECT_FALSE(seqan3::is_space('\0'));
}

TEST(char_predicate, is_space_msg)
{
    EXPECT_EQ((seqan3::is_space.message()), "(is_in_interval<'\t', '\r'> || is_char<' '>)"s);
}

TEST(char_predicate, is_punct)
{
    EXPECT_TRUE(seqan3::is_punct('!'));
    EXPECT_TRUE(seqan3::is_punct('"'));
    EXPECT_TRUE(seqan3::is_punct('.'));
    EXPECT_TRUE(seqan3::is_punct('/'));
    EXPECT_TRUE(seqan3::is_punct(':'));
    EXPECT_TRUE(seqan3::is_punct('@'));
    EXPECT_TRUE(seqan3::is_punct('['));
    EXPECT_TRUE(seqan3::is_punct('`'));
    EXPECT_TRUE(seqan3::is_punct('{'));
    EXPECT_TRUE(seqan3::is_punct('~'));
    EXPECT_FALSE(seqan3::is_punct(' '));
    EXPECT_FALSE(seqan3::is_punct('0'));
    EXPECT_FALSE(seqan3::is_punct('\0'));
}

TEST(char_predicate, is_punct_msg)
{
    EXPECT_EQ(
        (seqan3::is_punct.message()),
        "(((is_in_interval<'!', '/'> || is_in_interval<':', '@'>) || is_in_interval<'[', '`'>) || is_in_interval<'{', '~'>)"s);
}

TEST(char_predicate, is_alpha)
{
    EXPECT_FALSE(seqan3::is_alpha('\n'));
    EXPECT_FALSE(seqan3::is_alpha('\r'));
    EXPECT_FALSE(seqan3::is_alpha('\t'));
    EXPECT_FALSE(seqan3::is_alpha(' '));
    EXPECT_FALSE(seqan3::is_alpha('0'));
    EXPECT_TRUE(seqan3::is_alpha('a'));
    EXPECT_TRUE(seqan3::is_alpha('z'));
    EXPECT_TRUE(seqan3::is_alpha('Z'));
}

TEST(char_predicate, is_alpha_msg)
{
    EXPECT_EQ((seqan3::is_alpha.message()), "(is_in_interval<'A', 'Z'> || is_in_interval<'a', 'z'>)"s);
}

TEST(char_predicate, is_upper)
{
    EXPECT_FALSE(seqan3::is_upper('\n'));
    EXPECT_FALSE(seqan3::is_upper('\r'));
    EXPECT_FALSE(seqan3::is_upper('\t'));
    EXPECT_FALSE(seqan3::is_upper(' '));
    EXPECT_FALSE(seqan3::is_upper('0'));
    EXPECT_TRUE(seqan3::is_upper('A'));
    EXPECT_TRUE(seqan3::is_upper('Z'));
    EXPECT_FALSE(seqan3::is_upper('a'));
    EXPECT_FALSE(seqan3::is_upper('z'));
}

TEST(char_predicate, is_upper_msg)
{
    EXPECT_EQ((seqan3::is_upper.message()), "is_in_interval<'A', 'Z'>"s);
}

TEST(char_predicate, is_lower)
{
    EXPECT_FALSE(seqan3::is_lower('\n'));
    EXPECT_FALSE(seqan3::is_lower('\r'));
    EXPECT_FALSE(seqan3::is_lower('\t'));
    EXPECT_FALSE(seqan3::is_lower(' '));
    EXPECT_FALSE(seqan3::is_lower('0'));
    EXPECT_FALSE(seqan3::is_lower('A'));
    EXPECT_FALSE(seqan3::is_lower('Z'));
    EXPECT_TRUE(seqan3::is_lower('a'));
    EXPECT_TRUE(seqan3::is_lower('z'));
}

TEST(char_predicate, is_lower_msg)
{
    EXPECT_EQ((seqan3::is_lower.message()), "is_in_interval<'a', 'z'>"s);
}

TEST(char_predicate, is_digit)
{
    EXPECT_FALSE(seqan3::is_digit('\n'));
    EXPECT_FALSE(seqan3::is_digit('\r'));
    EXPECT_FALSE(seqan3::is_digit('\t'));
    EXPECT_FALSE(seqan3::is_digit(' '));
    EXPECT_TRUE(seqan3::is_digit('0'));
    EXPECT_TRUE(seqan3::is_digit('9'));
    EXPECT_FALSE(seqan3::is_digit('a'));
    EXPECT_FALSE(seqan3::is_digit('z'));
    EXPECT_FALSE(seqan3::is_digit('Z'));
}

TEST(char_predicate, is_digit_msg)
{
    EXPECT_EQ((seqan3::is_digit.message()), "is_in_interval<'0', '9'>"s);
}

TEST(char_predicate, is_xdigit)
{
    EXPECT_TRUE(seqan3::is_xdigit('0'));
    EXPECT_TRUE(seqan3::is_xdigit('9'));
    EXPECT_TRUE(seqan3::is_xdigit('a'));
    EXPECT_TRUE(seqan3::is_xdigit('f'));
    EXPECT_TRUE(seqan3::is_xdigit('A'));
    EXPECT_TRUE(seqan3::is_xdigit('F'));
    EXPECT_FALSE(seqan3::is_xdigit('g'));
    EXPECT_FALSE(seqan3::is_xdigit('z'));
    EXPECT_FALSE(seqan3::is_xdigit('G'));
    EXPECT_FALSE(seqan3::is_xdigit('Z'));
    EXPECT_FALSE(seqan3::is_xdigit('\n'));
    EXPECT_FALSE(seqan3::is_xdigit('\r'));
    EXPECT_FALSE(seqan3::is_xdigit('\t'));
    EXPECT_FALSE(seqan3::is_xdigit(' '));
}

TEST(char_predicate, is_xdigit_msg)
{
    EXPECT_EQ((seqan3::is_xdigit.message()),
              "((is_in_interval<'0', '9'> || is_in_interval<'A', 'F'>) || is_in_interval<'a', 'f'>)"s);
}

TEST(char_predicate, is_alnum)
{
    EXPECT_FALSE(seqan3::is_alnum('\n'));
    EXPECT_FALSE(seqan3::is_alnum('\r'));
    EXPECT_FALSE(seqan3::is_alnum('\t'));
    EXPECT_FALSE(seqan3::is_alnum(' '));
    EXPECT_TRUE(seqan3::is_alnum('0'));
    EXPECT_TRUE(seqan3::is_alnum('9'));
    EXPECT_TRUE(seqan3::is_alnum('a'));
    EXPECT_TRUE(seqan3::is_alnum('z'));
    EXPECT_TRUE(seqan3::is_alnum('Z'));
}

TEST(char_predicate, is_alnum_msg)
{
    EXPECT_EQ((seqan3::is_alnum.message()),
              "((is_in_interval<'0', '9'> || is_in_interval<'A', 'Z'>) || is_in_interval<'a', 'z'>)"s);
}

TEST(char_predicate, is_graph)
{
    EXPECT_FALSE(seqan3::is_graph('\n'));
    EXPECT_FALSE(seqan3::is_graph('\r'));
    EXPECT_FALSE(seqan3::is_graph('\t'));
    EXPECT_FALSE(seqan3::is_graph(' '));
    EXPECT_TRUE(seqan3::is_graph('0'));
    EXPECT_TRUE(seqan3::is_graph('9'));
    EXPECT_TRUE(seqan3::is_graph('a'));
    EXPECT_TRUE(seqan3::is_graph('z'));
    EXPECT_TRUE(seqan3::is_graph('Z'));
    EXPECT_TRUE(seqan3::is_graph('~'));
}

TEST(char_predicate, is_graph_msg)
{
    EXPECT_EQ((seqan3::is_graph.message()), "is_in_interval<'!', '~'>"s);
}

TEST(char_predicate, char_types)
{
    { // is_char
        char c1 = '\t';
        EXPECT_TRUE(seqan3::is_char<'\t'>(c1));
        char16_t c2 = '\t';
        EXPECT_TRUE(seqan3::is_char<'\t'>(c2));
        char32_t c3 = '\t';
        EXPECT_TRUE(seqan3::is_char<'\t'>(c3));
    }

    { // check value out of range.
        EXPECT_FALSE(seqan3::is_char<'\t'>(char16_t{256}));
    }

    { // is_in_interval
        char c1 = 'n';
        EXPECT_TRUE((seqan3::is_in_interval<'a', 'z'>(c1)));
        char16_t c2 = 'n';
        EXPECT_TRUE((seqan3::is_in_interval<'a', 'z'>(c2)));
        char32_t c3 = 'n';
        EXPECT_TRUE((seqan3::is_in_interval<'a', 'z'>(c3)));
    }

    { // check value out of range.
        EXPECT_FALSE((seqan3::is_in_interval<'a', 'z'>(char16_t{256})));
    }
}

// see issue https://github.com/seqan/seqan3/issues/1972
TEST(char_predicate, issue1972)
{
    EXPECT_TRUE(seqan3::char_is_valid_for<seqan3::gapped<seqan3::rna5>>('A'));  // valid seqan3::rna5 char
    EXPECT_TRUE(seqan3::char_is_valid_for<seqan3::gapped<seqan3::rna5>>('a'));  // valid seqan3::rna5 char
    EXPECT_TRUE(seqan3::char_is_valid_for<seqan3::gapped<seqan3::rna5>>('-'));  // valid seqan3::gap char
    EXPECT_FALSE(seqan3::char_is_valid_for<seqan3::gapped<seqan3::rna5>>('S')); // neither seqan3::rna5 nor seqan3::gap
}
