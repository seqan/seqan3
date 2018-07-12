// ==========================================================================
//                 SeqAn - The Library for Sequence Analysis
// ==========================================================================
//
// Copyright (c) 2006-2017, Knut Reinert, FU Berlin
// Copyright (c) 2016-2017, Knut Reinert & MPI Molekulare Genetik
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

#include <vector>
#include <array>
#include <list>
#include <forward_list>
#include <deque>
#include <string>

#include <seqan3/range/container/concept.hpp>
#include <seqan3/range/container/concatenated_sequences.hpp>

#include <sdsl/int_vector.hpp>

using namespace seqan3;

// if detail::sequence_container_concept_modified_by_const_iterator_bug<> is false
// test seqan3::concatenated_sequences<std::string>, otherwise
// test seqan3::concatenated_sequences<std::vector<char>>
using concatenated_sequences_string_t = seqan3::concatenated_sequences<
        std::conditional_t<
            detail::sequence_container_concept_modified_by_const_iterator_bug<>,
            std::vector<char>,
            std::string
      >>;

TEST(range_concept, forward_range_concept)
{
    EXPECT_TRUE((seqan3::forward_range_concept<std::array<char, 2>>));
    EXPECT_TRUE((seqan3::forward_range_concept<std::list<char>>));
    EXPECT_TRUE((seqan3::forward_range_concept<std::forward_list<char>>)); // `.size()` missing
    EXPECT_TRUE((seqan3::forward_range_concept<std::vector<char>>));
    EXPECT_TRUE((seqan3::forward_range_concept<std::deque<char>>));
    EXPECT_TRUE((seqan3::forward_range_concept<std::string>));

    EXPECT_TRUE((seqan3::forward_range_concept<concatenated_sequences_string_t>));
    EXPECT_TRUE((seqan3::forward_range_concept<seqan3::concatenated_sequences<std::vector<char>>>));
}

TEST(container_concept, container_concept)
{
    EXPECT_TRUE((seqan3::container_concept<std::array<char, 2>>));
    EXPECT_TRUE((seqan3::container_concept<std::list<char>>));
    EXPECT_FALSE((seqan3::container_concept<std::forward_list<char>>)); // `.size()` missing
    EXPECT_TRUE((seqan3::container_concept<std::vector<char>>));
    EXPECT_TRUE((seqan3::container_concept<std::deque<char>>));
    EXPECT_TRUE((seqan3::container_concept<std::string>));

    EXPECT_TRUE((seqan3::container_concept<concatenated_sequences_string_t>));
    EXPECT_TRUE((seqan3::container_concept<seqan3::concatenated_sequences<std::vector<char>>>));
}

template <typename string_t>
void container_concept_travis_bug_test()
{
    // non const version of container_concept_const_travis_bug_test
    using namespace std::string_literals;

    // example code from http://en.cppreference.com/w/cpp/string/basic_string/insert
    string_t s = "xmplr";
    string_t r = "";
    typename string_t::iterator it;

    // insert(size_type index, size_type count, char ch)
    r = s.insert(0, 1, 'E');
    EXPECT_EQ("Exmplr", s);

    // insert(size_type index, const char* s)
    r = s.insert(2, "e");
    EXPECT_EQ("Exemplr", s);

    // insert(size_type index, string const& str)
    r = s.insert(6, "a"s);
    EXPECT_EQ("Exemplar", s);

    // insert(size_type index, string const& str,
    //     size_type index_str, size_type count)
    r = s.insert(8, " is an example string."s, 0, 14);
    EXPECT_EQ("Exemplar is an example", s);

    // insert(const_iterator pos, char ch)
    it = s.insert(s.begin() + s.find_first_of('n') + 1, ':');
    EXPECT_EQ("Exemplar is an: example", s);

    // insert(const_iterator pos, size_type count, char ch)
    //TODO should return type::iterator on strings, remove if
    // sequence_container_concept_modified_by_const_iterator_bug is no issue anymore
    // it =
    s.insert(s.begin() + s.find_first_of(':') + 1, 2, '=');
    EXPECT_EQ("Exemplar is an:== example", s);

    // insert(const_iterator pos, InputIt first, InputIt last)
    {
        string_t seq = " string";
        //TODO should return type::iterator on strings, remove if
        // sequence_container_concept_modified_by_const_iterator_bug is no issue anymore
        // it =
        s.insert(s.begin() + s.find_last_of('e') + 1,
            std::begin(seq), std::end(seq));
        EXPECT_EQ("Exemplar is an:== example string", s);
    }

    // insert(const_iterator pos, std::initializer_list<char>)
    //TODO should return type::iterator on strings (https://gcc.gnu.org/bugzilla/show_bug.cgi?id=83328)
    // it =
    s.insert(s.begin() + s.find_first_of('g') + 1, { '.' });
    EXPECT_EQ("Exemplar is an:== example string.", s);
}

template <typename string_t>
void container_concept_const_travis_bug_test()
{
    // travis failed on this statement
    // concept bool sequence_container_concept = requires (type val, type val2)
    //              ^~~~~~~~~~~~~~~~
    // /include/seqan3/range/container/concept.hpp:113:14: note:     with ‘std::basic_string<char> val’
    // /include/seqan3/range/container/concept.hpp:113:14: note:     with ‘std::basic_string<char> val2’
    // [...]
    // /include/seqan3/range/container/concept.hpp:113:14: note: the required expression ‘val.erase(val.cbegin(), val.cend())’ would be ill-formed

    using namespace std::string_literals;

    static_assert(detail::sequence_container_concept_modified_by_const_iterator<string_t>);
    static_assert(!detail::sequence_container_concept_modified_by_const_iterator_bug<string_t> && std::is_same_v<string_t, std::string>);
    static_assert(seqan3::container_concept<string_t>);

    // example code from http://en.cppreference.com/w/cpp/string/basic_string/insert
    string_t s = "xmplr";

    // insert(size_type index, size_type count, char ch)
    s.insert(0, 1, 'E');
    EXPECT_EQ("Exmplr", s);

    // insert(size_type index, const char* s)
    s.insert(2, "e");
    EXPECT_EQ("Exemplr", s);

    // insert(size_type index, string const& str)
    s.insert(6, "a"s);
    EXPECT_EQ("Exemplar", s);

    // insert(size_type index, string const& str,
    //     size_type index_str, size_type count)
    s.insert(8, " is an example string."s, 0, 14);
    EXPECT_EQ("Exemplar is an example", s);

    // insert(const_iterator pos, char ch)
    s.insert(s.cbegin() + s.find_first_of('n') + 1, ':');
    EXPECT_EQ("Exemplar is an: example", s);

    // insert(const_iterator pos, size_type count, char ch)
    s.insert(s.cbegin() + s.find_first_of(':') + 1, 2, '=');
    EXPECT_EQ("Exemplar is an:== example", s);

    // insert(const_iterator pos, InputIt first, InputIt last)
    {
        string_t seq = " string";
        s.insert(s.cbegin() + s.find_last_of('e') + 1,
            std::begin(seq), std::end(seq));
        EXPECT_EQ("Exemplar is an:== example string", s);
    }

    // insert(const_iterator pos, std::initializer_list<char>)
    s.insert(s.cbegin() + s.find_first_of('g') + 1, { '.' });
    EXPECT_EQ("Exemplar is an:== example string.", s);
}

TEST(container_concept, container_concept_travis_bug)
{
    container_concept_travis_bug_test<std::string>();

    if constexpr(!detail::sequence_container_concept_modified_by_const_iterator_bug<>)
        container_concept_const_travis_bug_test<std::string>();
}

TEST(container_concept, sequence_container_concept)
{
    EXPECT_FALSE((seqan3::sequence_container_concept<std::array<char, 2>>));
    EXPECT_TRUE((seqan3::sequence_container_concept<std::list<char>>));
    EXPECT_FALSE((seqan3::sequence_container_concept<std::forward_list<char>>));
    EXPECT_TRUE((seqan3::sequence_container_concept<std::vector<char>>));
    EXPECT_TRUE((seqan3::sequence_container_concept<std::deque<char>>));
    EXPECT_TRUE((seqan3::sequence_container_concept<std::string>));

    EXPECT_TRUE((seqan3::sequence_container_concept<concatenated_sequences_string_t>));
    EXPECT_TRUE((seqan3::sequence_container_concept<seqan3::concatenated_sequences<std::vector<char>>>));
}

TEST(container_concept, sequence_container_concept_modified_by_const_iterator)
{
    EXPECT_FALSE((seqan3::detail::sequence_container_concept_modified_by_const_iterator<std::array<char, 2>>));
    EXPECT_TRUE((seqan3::detail::sequence_container_concept_modified_by_const_iterator<std::list<char>>));
    EXPECT_FALSE((seqan3::detail::sequence_container_concept_modified_by_const_iterator<std::forward_list<char>>));
    EXPECT_TRUE((seqan3::detail::sequence_container_concept_modified_by_const_iterator<std::vector<char>>));
    EXPECT_TRUE((seqan3::detail::sequence_container_concept_modified_by_const_iterator<std::deque<char>>));
    if constexpr(!detail::sequence_container_concept_modified_by_const_iterator_bug<>)
    {
        EXPECT_TRUE((seqan3::detail::sequence_container_concept_modified_by_const_iterator<std::string>));
    }

    EXPECT_TRUE((seqan3::detail::sequence_container_concept_modified_by_const_iterator<concatenated_sequences_string_t>));
    EXPECT_TRUE((seqan3::detail::sequence_container_concept_modified_by_const_iterator<seqan3::concatenated_sequences<std::vector<char>>>));
}

TEST(container_concept, random_access_container_concept)
{
    EXPECT_FALSE((seqan3::random_access_container_concept<std::array<char, 2>>));
    EXPECT_FALSE((seqan3::random_access_container_concept<std::list<char>>));
    EXPECT_FALSE((seqan3::random_access_container_concept<std::forward_list<char>>));
    EXPECT_TRUE((seqan3::random_access_container_concept<std::vector<char>>));
    EXPECT_TRUE((seqan3::random_access_container_concept<std::deque<char>>));
    EXPECT_TRUE((seqan3::random_access_container_concept<std::string>));

    EXPECT_TRUE((seqan3::random_access_container_concept<concatenated_sequences_string_t>));
    EXPECT_TRUE((seqan3::random_access_container_concept<seqan3::concatenated_sequences<std::vector<char>>>));
}

TEST(container_concept, reservable_container_concept)
{
    EXPECT_FALSE((seqan3::reservable_container_concept<std::array<char, 2>>));
    EXPECT_FALSE((seqan3::reservable_container_concept<std::list<char>>));
    EXPECT_FALSE((seqan3::reservable_container_concept<std::forward_list<char>>));
    EXPECT_TRUE((seqan3::reservable_container_concept<std::vector<char>>));
    EXPECT_FALSE((seqan3::reservable_container_concept<std::deque<char>>));
    EXPECT_TRUE((seqan3::reservable_container_concept<std::string>));

    EXPECT_TRUE((seqan3::reservable_container_concept<concatenated_sequences_string_t>));
    EXPECT_TRUE((seqan3::reservable_container_concept<seqan3::concatenated_sequences<std::vector<char>>>));

    EXPECT_TRUE((seqan3::reservable_container_concept<sdsl::bit_vector>));
    EXPECT_TRUE((seqan3::reservable_container_concept<sdsl::int_vector<>>));
    EXPECT_TRUE((seqan3::reservable_container_concept<sdsl::int_vector<13>>));
    EXPECT_TRUE((seqan3::reservable_container_concept<sdsl::int_vector<64>>));
}

/* Check the SDSL containers */
//TODO

/* Check our containers */
//TODO
