// ==========================================================================
//                 SeqAn - The Library for Sequence Analysis
// ==========================================================================
//
// Copyright (c) 2006-2018, Knut Reinert, FU Berlin
// Copyright (c) 2016-2018, Knut Reinert & MPI Molekulare Genetik
// All rights reserved.
//
// Redistribution and use in source and binary forms, with or without
// modification, are permitted provided that the following conditions are met:
//
//     * Redistributions of source code must retain the above copyright
//       notice, this list of conditions and the following disclaimer.
//     * Redistributions in binary form must reproduce the above copyright
//       notice, this list of conditions and the following disclaimer in the
//       documentation and/or other materials provided with the distribution.
//     * Neither the name of Knut Reinert or the FU Berlin nor the names of
//       its contributors may be used to endorse or promote products derived
//       from this software without specific prior written permission.
//
// THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
// AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
// IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
// ARE DISCLAIMED. IN NO EVENT SHALL KNUT REINERT OR THE FU BERLIN BE LIABLE
// FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL
// DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR
// SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER
// CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT
// LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY
// OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH
// DAMAGE.
//
// ==========================================================================

#include <gtest/gtest.h>

#include <string>

#include <seqan3/alphabet/all.hpp>
#include <seqan3/io/stream/parse_condition.hpp>

using namespace std::literals;
using namespace seqan3;

template <char char_v>
struct foo : seqan3::detail::parse_condition_base<foo<char_v>>
{
    inline static constexpr constexpr_string msg{constexpr_string{"foo_"} + constexpr_string{char_v}};

    using base_t = seqan3::detail::parse_condition_base<foo<char_v>>;

    static auto constexpr data = [&] () { typename base_t::data_t d{}; d[char_v] = true; return d; }();
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

TEST(parse_condition, parse_condition)
{
    foo<'a'> p{};
    EXPECT_TRUE(p('a'));
    EXPECT_FALSE(p('f'));
}

TEST(parse_condition, parse_condition_msg)
{
    EXPECT_EQ(foo<'o'>::msg.string(), "foo_o"s);
}

TEST(parse_condition, parse_condition_concept)
{
    using namespace seqan3;

    EXPECT_TRUE(detail::parse_condition_concept<decltype(is_in_alphabet<dna4>)>);
    EXPECT_TRUE(detail::parse_condition_concept<decltype(is_char<to_char(aa27::A)>)>);
    EXPECT_TRUE((detail::parse_condition_concept<decltype(is_in_interval<'a','z'>)>));
    EXPECT_TRUE(detail::parse_condition_concept<decltype(is_space)>);
    EXPECT_TRUE(detail::parse_condition_concept<decltype(is_blank)>);
    EXPECT_TRUE(detail::parse_condition_concept<decltype(is_graph)>);
    EXPECT_TRUE(detail::parse_condition_concept<decltype(is_alpha)>);
    EXPECT_TRUE(detail::parse_condition_concept<decltype(is_digit)>);
    EXPECT_TRUE(detail::parse_condition_concept<decltype(is_alnum)>);
    // TODO(rrahn): Direct call in concept caused ICE.
    auto val = ((!is_space || is_alpha) || is_digit);
    EXPECT_TRUE(detail::parse_condition_concept<decltype(val)>);
    EXPECT_TRUE(detail::parse_condition_concept<foo<' '>>);
    EXPECT_FALSE(detail::parse_condition_concept<bar>);
    EXPECT_FALSE(detail::parse_condition_concept<int>);
}

TEST(parse_condition, parse_condition_combiner)
{
    using namespace seqan3::detail;

    using cond_t = detail::parse_condition_combiner<foo<'a'>, foo<'A'>, foo<'0'>>;
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

TEST(parse_condition, parse_condition_combiner_msg)
{
    using namespace seqan3::detail;
    using or_t = detail::parse_condition_combiner<foo<'a'>, foo<'A'>, foo<'0'>>;
    EXPECT_EQ(or_t::msg.string(),   "(foo_a || foo_A || foo_0)"s);
}

TEST(parse_condition, is_not)
{
    using namespace seqan3::detail;
    using cond_t = detail::parse_condition_negator<foo<'a'>>;
    EXPECT_FALSE(cond_t{}('a'));
    EXPECT_TRUE(cond_t{}('A'));
    EXPECT_TRUE(cond_t{}('0'));

    auto constexpr p = !foo_v<'a'>;
    EXPECT_FALSE(p('a'));
    EXPECT_TRUE(p('A'));
    EXPECT_TRUE(p('0'));
}

TEST(parse_condition, is_not_msg)
{
    using namespace seqan3::detail;
    using fn = decltype(!is_alpha);
    EXPECT_EQ(fn::msg.string(), "!(is_in_interval<'A', 'Z'> || is_in_interval<'a', 'z'>)"s);
}

TEST(parse_condition, is_in_interval)
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

TEST(parse_condition, is_in_interval_msg)
{
    using namespace seqan3;
    EXPECT_EQ((detail::is_in_interval_type<'a', 'z'>::msg.string()), "is_in_interval<'a', 'z'>"s);
}

TEST(parse_condition, is_in_alphabet)
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

TEST(parse_condition, is_in_alphabet_msg)
{
    using namespace seqan3;
    EXPECT_EQ((detail::is_in_alphabet_type<dna4>::msg.string()), "is_in_alphabet<seqan3::dna4>"s);
}

TEST(parse_condition, is_char)
{
    using namespace seqan3;
    {
        auto constexpr cond = is_char<'A'>;
        EXPECT_TRUE(cond('A'));
        EXPECT_FALSE(cond('x'));
    }

    {
        auto constexpr cond = is_char<to_char(aa27::A)>;
        EXPECT_TRUE(cond('A'));
        EXPECT_FALSE(cond('z'));
    }
}

TEST(parse_condition, is_char_msg)
{
    using namespace seqan3;
    EXPECT_EQ((is_char<to_char('A'_dna4)>.msg.string()), "is_char<'A'>"s);
    EXPECT_EQ((is_char<'\t'>.msg.string()), "is_char<'\t'>"s);
}

TEST(parse_condition, is_cntrl)
{
    using namespace seqan3;
    EXPECT_TRUE(is_cntrl('\0'));
    EXPECT_TRUE(is_cntrl(static_cast<char>(31)));
    EXPECT_TRUE(is_cntrl(static_cast<char>(127)));
    EXPECT_TRUE(is_cntrl('\t'));
    EXPECT_FALSE(is_cntrl('A'));
}

TEST(parse_condition, is_print)
{
    using namespace seqan3;
    EXPECT_FALSE(is_print('\0'));
    EXPECT_FALSE(is_print(static_cast<char>(31)));
    EXPECT_FALSE(is_print(static_cast<char>(127)));
    EXPECT_TRUE(is_print(' '));
    EXPECT_TRUE(is_print('A'));
    EXPECT_TRUE(is_print('~'));
}

TEST(parse_condition, is_print_msg)
{
    using namespace seqan3;
    EXPECT_EQ((is_print.message()), "is_in_interval<' ', '~'>"s);
}


TEST(parse_condition, is_blank)
{
    using namespace seqan3;
    EXPECT_TRUE(is_blank(' '));
    EXPECT_TRUE(is_blank('\t'));
    EXPECT_FALSE(is_blank('A'));
    EXPECT_FALSE(is_blank('\n'));
}

TEST(parse_condition, is_blank_msg)
{
    using namespace seqan3;
    EXPECT_EQ((is_blank.message()), "(is_char<'\t'> || is_char<' '>)"s);
}

TEST(parse_condition, is_space)
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

TEST(parse_condition, is_space_msg)
{
    using namespace seqan3;
    EXPECT_EQ((is_space.message()),
              "(is_in_interval<'\t', '\r'> || is_char<' '>)"s);
}

TEST(parse_condition, is_punct)
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

TEST(parse_condition, is_punct_msg)
{
    using namespace seqan3;
    EXPECT_EQ((is_punct.message()),
              "(((is_in_interval<'!', '/'> || is_in_interval<':', '@'>) || is_in_interval<'[', '`'>) || is_in_interval<'{', '~'>)"s);
}

TEST(parse_condition, is_alpha)
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

TEST(parse_condition, is_alpha_msg)
{
    using namespace seqan3;
    EXPECT_EQ((is_alpha.message()),
              "(is_in_interval<'A', 'Z'> || is_in_interval<'a', 'z'>)"s);
}

TEST(parse_condition, is_upper)
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

TEST(parse_condition, is_upper_msg)
{
    using namespace seqan3;
    EXPECT_EQ((is_upper.message()), "is_in_interval<'A', 'Z'>"s);
}

TEST(parse_condition, is_lower)
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

TEST(parse_condition, is_lower_msg)
{
    using namespace seqan3;
    EXPECT_EQ((is_lower.message()), "is_in_interval<'a', 'z'>"s);
}

TEST(parse_condition, is_digit)
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

TEST(parse_condition, is_digit_msg)
{
    using namespace seqan3;
    EXPECT_EQ((is_digit.message()), "is_in_interval<'0', '9'>"s);
}

TEST(parse_condition, is_xdigit)
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

TEST(parse_condition, is_xdigit_msg)
{
    using namespace seqan3;
    EXPECT_EQ((is_xdigit.message()), "((is_in_interval<'0', '9'> || is_in_interval<'A', 'F'>) || is_in_interval<'a', 'f'>)"s);
}

TEST(parse_condition, is_alnum)
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

TEST(parse_condition, is_alnum_msg)
{
    using namespace seqan3;
    EXPECT_EQ((is_alnum.message()), "((is_in_interval<'0', '9'> || is_in_interval<'A', 'Z'>) || is_in_interval<'a', 'z'>)"s);
}

TEST(parse_condition, is_graph)
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

TEST(parse_condition, is_graph_msg)
{
    using namespace seqan3;
    EXPECT_EQ((is_graph.message()), "is_in_interval<'!', '~'>"s);
}

TEST(parse_condition, char_types)
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

TEST(parse_condition, parse_asserter)
{
    using namespace seqan3;

    parse_asserter asserter{is_alnum};
    EXPECT_THROW(asserter('\t'), parse_error);
    EXPECT_NO_THROW(asserter('a'));

    try
    {
        asserter('\t');
    } catch(parse_error & e)
    {
        EXPECT_EQ(e.what(), "Parsed value <'\\t'> which does not fulfill the following condition:"
                            " ((is_in_interval<'0', '9'> || is_in_interval<'A', 'Z'>) || is_in_interval<'a', 'z'>)"s);
    }
}
