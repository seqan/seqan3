// -----------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2019, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2019, Knut Reinert & MPI für molekulare Genetik
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at: https://github.com/seqan/seqan3/blob/master/LICENSE.md
// -----------------------------------------------------------------------------------------------------

#include <gtest/gtest.h>

#include <string>

#include <seqan3/alphabet/all.hpp>
#include <seqan3/core/char_operations/predicate.hpp>

using namespace std::literals;
using namespace seqan3;

template <char char_v>
struct foo : seqan3::detail::char_predicate_base<foo<char_v>>
{
    inline static constexpr small_string msg{small_string{"foo_"} + small_string{char_v}};

    using base_t = seqan3::detail::char_predicate_base<foo<char_v>>;

    static auto constexpr data = [] () { typename base_t::data_t d{}; d[char_v] = true; return d; }();
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

TEST(char_predicate, char_predicate)
{
    foo<'a'> p{};
    EXPECT_TRUE(p('a'));
    EXPECT_FALSE(p('f'));
}

TEST(char_predicate, char_predicate_msg)
{
    EXPECT_EQ(foo<'o'>::msg.str(), "foo_o"s);
}

TEST(char_predicate, CharPredicate)
{
    using namespace seqan3;

    EXPECT_TRUE(detail::CharPredicate<decltype(is_in_alphabet<dna4>)>);
    EXPECT_TRUE(detail::CharPredicate<decltype(is_char<to_char('A'_aa27)>)>);
    EXPECT_TRUE((detail::CharPredicate<decltype(is_in_interval<'a','z'>)>));
    EXPECT_TRUE(detail::CharPredicate<decltype(is_space)>);
    EXPECT_TRUE(detail::CharPredicate<decltype(is_blank)>);
    EXPECT_TRUE(detail::CharPredicate<decltype(is_graph)>);
    EXPECT_TRUE(detail::CharPredicate<decltype(is_alpha)>);
    EXPECT_TRUE(detail::CharPredicate<decltype(is_digit)>);
    EXPECT_TRUE(detail::CharPredicate<decltype(is_alnum)>);
    // TODO(rrahn): Direct call in concept caused ICE.
    auto val = ((!is_space || is_alpha) || is_digit);
    EXPECT_TRUE(detail::CharPredicate<decltype(val)>);
    EXPECT_TRUE(detail::CharPredicate<foo<' '>>);
    EXPECT_FALSE(detail::CharPredicate<bar>);
    EXPECT_FALSE(detail::CharPredicate<int>);
}

TEST(char_predicate, char_predicate_combiner)
{
    using namespace seqan3::detail;

    using cond_t = detail::char_predicate_combiner<foo<'a'>, foo<'A'>, foo<'0'>>;
    EXPECT_TRUE(cond_t{}('a'));
    EXPECT_TRUE(cond_t{}('A'));
    EXPECT_TRUE(cond_t{}('0'));
    EXPECT_FALSE(cond_t{}('z'));
    EXPECT_FALSE(cond_t{}('!'));
    EXPECT_FALSE(cond_t{}('1'));

    auto constexpr p = foo_v<'a'> || foo_v<'A'> || foo_v<'0'>;
    EXPECT_TRUE(p('a'));
    EXPECT_TRUE(p('A'));
    EXPECT_TRUE(p('0'));
    EXPECT_FALSE(p('z'));
    EXPECT_FALSE(p('!'));
    EXPECT_FALSE(p('1'));
}

TEST(char_predicate, char_predicate_combiner_msg)
{
    using namespace seqan3::detail;
    using or_t = detail::char_predicate_combiner<foo<'a'>, foo<'A'>, foo<'0'>>;
    EXPECT_EQ(or_t::msg.str(),   "(foo_a || foo_A || foo_0)"s);
}

TEST(char_predicate, is_not)
{
    using namespace seqan3::detail;
    using cond_t = detail::char_predicate_negator<foo<'a'>>;
    EXPECT_FALSE(cond_t{}('a'));
    EXPECT_TRUE(cond_t{}('A'));
    EXPECT_TRUE(cond_t{}('0'));

    auto constexpr p = !foo_v<'a'>;
    EXPECT_FALSE(p('a'));
    EXPECT_TRUE(p('A'));
    EXPECT_TRUE(p('0'));
}

TEST(char_predicate, is_not_msg)
{
    using namespace seqan3::detail;
    using fn = decltype(!is_alpha);
    EXPECT_EQ(fn::msg.str(), "!(is_in_interval<'A', 'Z'> || is_in_interval<'a', 'z'>)"s);
}

TEST(char_predicate, is_in_interval)
{
    using namespace seqan3;
    auto constexpr cond = is_in_interval<'a', 'z'>;
    EXPECT_TRUE(cond('a'));
    EXPECT_TRUE(cond('k'));
    EXPECT_TRUE(cond('z'));
    EXPECT_FALSE(cond('A'));
    EXPECT_FALSE(cond('0'));
    EXPECT_FALSE(cond('!'));
}

TEST(char_predicate, is_in_interval_msg)
{
    using namespace seqan3;
    EXPECT_EQ((detail::is_in_interval_type<'a', 'z'>::msg.str()), "is_in_interval<'a', 'z'>"s);
}

TEST(char_predicate, is_in_alphabet)
{
    using namespace seqan3;
    {
        auto constexpr cond = is_in_alphabet<dna4>;
        EXPECT_TRUE(cond('a'));
        EXPECT_TRUE(cond('A'));
        EXPECT_TRUE(cond('c'));
        EXPECT_TRUE(cond('C'));
        EXPECT_TRUE(cond('g'));
        EXPECT_TRUE(cond('G'));
        EXPECT_TRUE(cond('t'));
        EXPECT_TRUE(cond('T'));
        EXPECT_FALSE(cond('N'));
        EXPECT_FALSE(cond('n'));
        EXPECT_FALSE(cond('!'));
        EXPECT_FALSE(cond('0'));
    }

    {
        auto constexpr cond = is_in_alphabet<aa27>;
        EXPECT_TRUE(cond('a'));
        EXPECT_TRUE(cond('A'));
        EXPECT_TRUE(cond('z'));
        EXPECT_TRUE(cond('Z'));
        EXPECT_TRUE(cond('*'));
        EXPECT_FALSE(cond('!'));
        EXPECT_FALSE(cond('0'));
    }
}

TEST(char_predicate, is_in_alphabet_msg)
{
    using namespace seqan3;
    EXPECT_EQ((detail::is_in_alphabet_type<dna4>::msg.str()), "is_in_alphabet<seqan3::dna4>"s);
}

TEST(char_predicate, is_char)
{
    using namespace seqan3;
    {
        auto constexpr cond = is_char<'A'>;
        EXPECT_TRUE(cond('A'));
        EXPECT_FALSE(cond('x'));
    }

    {
        auto constexpr cond = is_char<to_char('A'_aa27)>;
        EXPECT_TRUE(cond('A'));
        EXPECT_FALSE(cond('z'));
    }
}

TEST(char_predicate, is_char_msg)
{
    using namespace seqan3;
    EXPECT_EQ((is_char<to_char('A'_dna4)>.msg.str()), "is_char<'A'>"s);
    EXPECT_EQ((is_char<'\t'>.msg.str()), "is_char<'\t'>"s);
}

TEST(char_predicate, is_cntrl)
{
    using namespace seqan3;
    EXPECT_TRUE(is_cntrl('\0'));
    EXPECT_TRUE(is_cntrl(static_cast<char>(31)));
    EXPECT_TRUE(is_cntrl(static_cast<char>(127)));
    EXPECT_TRUE(is_cntrl('\t'));
    EXPECT_FALSE(is_cntrl('A'));
}

TEST(char_predicate, is_print)
{
    using namespace seqan3;
    EXPECT_FALSE(is_print('\0'));
    EXPECT_FALSE(is_print(static_cast<char>(31)));
    EXPECT_FALSE(is_print(static_cast<char>(127)));
    EXPECT_TRUE(is_print(' '));
    EXPECT_TRUE(is_print('A'));
    EXPECT_TRUE(is_print('~'));
}

TEST(char_predicate, is_print_msg)
{
    using namespace seqan3;
    EXPECT_EQ((is_print.message()), "is_in_interval<' ', '~'>"s);
}


TEST(char_predicate, is_blank)
{
    using namespace seqan3;
    EXPECT_TRUE(is_blank(' '));
    EXPECT_TRUE(is_blank('\t'));
    EXPECT_FALSE(is_blank('A'));
    EXPECT_FALSE(is_blank('\n'));
}

TEST(char_predicate, is_blank_msg)
{
    using namespace seqan3;
    EXPECT_EQ((is_blank.message()), "(is_char<'\t'> || is_char<' '>)"s);
}

TEST(char_predicate, is_space)
{
    using namespace seqan3;
    EXPECT_TRUE(is_space('\n'));
    EXPECT_TRUE(is_space('\r'));
    EXPECT_TRUE(is_space('\f'));
    EXPECT_TRUE(is_space('\v'));
    EXPECT_TRUE(is_space('\t'));
    EXPECT_TRUE(is_space(' '));
    EXPECT_FALSE(is_space('0'));
    EXPECT_FALSE(is_space('\0'));
}

TEST(char_predicate, is_space_msg)
{
    using namespace seqan3;
    EXPECT_EQ((is_space.message()),
              "(is_in_interval<'\t', '\r'> || is_char<' '>)"s);
}

TEST(char_predicate, is_punct)
{
    using namespace seqan3;
    EXPECT_TRUE(is_punct('!'));
    EXPECT_TRUE(is_punct('"'));
    EXPECT_TRUE(is_punct('.'));
    EXPECT_TRUE(is_punct('/'));
    EXPECT_TRUE(is_punct(':'));
    EXPECT_TRUE(is_punct('@'));
    EXPECT_TRUE(is_punct('['));
    EXPECT_TRUE(is_punct('`'));
    EXPECT_TRUE(is_punct('{'));
    EXPECT_TRUE(is_punct('~'));
    EXPECT_FALSE(is_punct(' '));
    EXPECT_FALSE(is_punct('0'));
    EXPECT_FALSE(is_punct('\0'));
}

TEST(char_predicate, is_punct_msg)
{
    using namespace seqan3;
    EXPECT_EQ((is_punct.message()),
              "(((is_in_interval<'!', '/'> || is_in_interval<':', '@'>) || is_in_interval<'[', '`'>) || is_in_interval<'{', '~'>)"s);
}

TEST(char_predicate, is_alpha)
{
    using namespace seqan3;
    EXPECT_FALSE(is_alpha('\n'));
    EXPECT_FALSE(is_alpha('\r'));
    EXPECT_FALSE(is_alpha('\t'));
    EXPECT_FALSE(is_alpha(' '));
    EXPECT_FALSE(is_alpha('0'));
    EXPECT_TRUE(is_alpha('a'));
    EXPECT_TRUE(is_alpha('z'));
    EXPECT_TRUE(is_alpha('Z'));
}

TEST(char_predicate, is_alpha_msg)
{
    using namespace seqan3;
    EXPECT_EQ((is_alpha.message()),
              "(is_in_interval<'A', 'Z'> || is_in_interval<'a', 'z'>)"s);
}

TEST(char_predicate, is_upper)
{
    using namespace seqan3;
    EXPECT_FALSE(is_upper('\n'));
    EXPECT_FALSE(is_upper('\r'));
    EXPECT_FALSE(is_upper('\t'));
    EXPECT_FALSE(is_upper(' '));
    EXPECT_FALSE(is_upper('0'));
    EXPECT_TRUE(is_upper('A'));
    EXPECT_TRUE(is_upper('Z'));
    EXPECT_FALSE(is_upper('a'));
    EXPECT_FALSE(is_upper('z'));
}

TEST(char_predicate, is_upper_msg)
{
    using namespace seqan3;
    EXPECT_EQ((is_upper.message()), "is_in_interval<'A', 'Z'>"s);
}

TEST(char_predicate, is_lower)
{
    using namespace seqan3;
    EXPECT_FALSE(is_lower('\n'));
    EXPECT_FALSE(is_lower('\r'));
    EXPECT_FALSE(is_lower('\t'));
    EXPECT_FALSE(is_lower(' '));
    EXPECT_FALSE(is_lower('0'));
    EXPECT_FALSE(is_lower('A'));
    EXPECT_FALSE(is_lower('Z'));
    EXPECT_TRUE(is_lower('a'));
    EXPECT_TRUE(is_lower('z'));
}

TEST(char_predicate, is_lower_msg)
{
    using namespace seqan3;
    EXPECT_EQ((is_lower.message()), "is_in_interval<'a', 'z'>"s);
}

TEST(char_predicate, is_digit)
{
    using namespace seqan3;
    EXPECT_FALSE(is_digit('\n'));
    EXPECT_FALSE(is_digit('\r'));
    EXPECT_FALSE(is_digit('\t'));
    EXPECT_FALSE(is_digit(' '));
    EXPECT_TRUE(is_digit('0'));
    EXPECT_TRUE(is_digit('9'));
    EXPECT_FALSE(is_digit('a'));
    EXPECT_FALSE(is_digit('z'));
    EXPECT_FALSE(is_digit('Z'));
}

TEST(char_predicate, is_digit_msg)
{
    using namespace seqan3;
    EXPECT_EQ((is_digit.message()), "is_in_interval<'0', '9'>"s);
}

TEST(char_predicate, is_xdigit)
{
    using namespace seqan3;

    EXPECT_TRUE(is_xdigit('0'));
    EXPECT_TRUE(is_xdigit('9'));
    EXPECT_TRUE(is_xdigit('a'));
    EXPECT_TRUE(is_xdigit('f'));
    EXPECT_TRUE(is_xdigit('A'));
    EXPECT_TRUE(is_xdigit('F'));
    EXPECT_FALSE(is_xdigit('g'));
    EXPECT_FALSE(is_xdigit('z'));
    EXPECT_FALSE(is_xdigit('G'));
    EXPECT_FALSE(is_xdigit('Z'));
    EXPECT_FALSE(is_xdigit('\n'));
    EXPECT_FALSE(is_xdigit('\r'));
    EXPECT_FALSE(is_xdigit('\t'));
    EXPECT_FALSE(is_xdigit(' '));
}

TEST(char_predicate, is_xdigit_msg)
{
    using namespace seqan3;
    EXPECT_EQ((is_xdigit.message()), "((is_in_interval<'0', '9'> || is_in_interval<'A', 'F'>) || is_in_interval<'a', 'f'>)"s);
}

TEST(char_predicate, is_alnum)
{
    using namespace seqan3;
    EXPECT_FALSE(is_alnum('\n'));
    EXPECT_FALSE(is_alnum('\r'));
    EXPECT_FALSE(is_alnum('\t'));
    EXPECT_FALSE(is_alnum(' '));
    EXPECT_TRUE(is_alnum('0'));
    EXPECT_TRUE(is_alnum('9'));
    EXPECT_TRUE(is_alnum('a'));
    EXPECT_TRUE(is_alnum('z'));
    EXPECT_TRUE(is_alnum('Z'));
}

TEST(char_predicate, is_alnum_msg)
{
    using namespace seqan3;
    EXPECT_EQ((is_alnum.message()), "((is_in_interval<'0', '9'> || is_in_interval<'A', 'Z'>) || is_in_interval<'a', 'z'>)"s);
}

TEST(char_predicate, is_graph)
{
    using namespace seqan3;
    EXPECT_FALSE(is_graph('\n'));
    EXPECT_FALSE(is_graph('\r'));
    EXPECT_FALSE(is_graph('\t'));
    EXPECT_FALSE(is_graph(' '));
    EXPECT_TRUE(is_graph('0'));
    EXPECT_TRUE(is_graph('9'));
    EXPECT_TRUE(is_graph('a'));
    EXPECT_TRUE(is_graph('z'));
    EXPECT_TRUE(is_graph('Z'));
    EXPECT_TRUE(is_graph('~'));
}

TEST(char_predicate, is_graph_msg)
{
    using namespace seqan3;
    EXPECT_EQ((is_graph.message()), "is_in_interval<'!', '~'>"s);
}

TEST(char_predicate, char_types)
{
    using namespace seqan3;

    {  // is_char
        char c1 = '\t';
        EXPECT_TRUE(is_char<'\t'>(c1));
        char16_t c2 = '\t';
        EXPECT_TRUE(is_char<'\t'>(c2));
        char32_t c3 = '\t';
        EXPECT_TRUE(is_char<'\t'>(c3));
    }

    {  // check value out of range.
        EXPECT_FALSE(is_char<'\t'>(char16_t{256}));
    }

    {  // is_in_interval
        char c1 = 'n';
        EXPECT_TRUE((is_in_interval<'a', 'z'>(c1)));
        char16_t c2 = 'n';
        EXPECT_TRUE((is_in_interval<'a', 'z'>(c2)));
        char32_t c3 = 'n';
        EXPECT_TRUE((is_in_interval<'a', 'z'>(c3)));
    }

    {  // check value out of range.
        EXPECT_FALSE((is_in_interval<'a', 'z'>(char16_t{256})));
    }

    {  // is_in_alphabet
        char c1 = 'N';
        EXPECT_TRUE(is_in_alphabet<seqan3::dna5>(c1));
        char16_t c2 = 'N';
        EXPECT_TRUE(is_in_alphabet<seqan3::dna5>(c2));
        char32_t c3 = 'N';
        EXPECT_TRUE(is_in_alphabet<seqan3::dna5>(c3));
    }

    {  // check value out of range
        EXPECT_FALSE(is_in_alphabet<seqan3::dna5>(char16_t{256}));
    }
}
