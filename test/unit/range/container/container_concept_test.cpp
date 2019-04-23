// -----------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2019, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2019, Knut Reinert & MPI für molekulare Genetik
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at: https://github.com/seqan/seqan3/blob/master/LICENSE.md
// -----------------------------------------------------------------------------------------------------

#include <gtest/gtest.h>

#include <vector>
#include <array>
#include <list>
#include <forward_list>
#include <deque>
#include <string>

#include <sdsl/int_vector.hpp>

#include <seqan3/range/container/concatenated_sequences.hpp>
#include <seqan3/range/container/concept.hpp>

using namespace seqan3;

// if detail::SequenceContainer_modified_by_const_iterator_bug<> is false
// test seqan3::concatenated_sequences<std::string>, otherwise
// test seqan3::concatenated_sequences<std::vector<char>>
using concatenated_sequences_string_t = seqan3::concatenated_sequences<
        std::conditional_t<
            detail::SequenceContainer_modified_by_const_iterator_bug<>,
            std::vector<char>,
            std::string
      >>;

TEST(container_concept, ForwardRange)
{
    EXPECT_TRUE((std::ranges::ForwardRange<std::array<char, 2>>));
    EXPECT_TRUE((std::ranges::ForwardRange<std::list<char>>));
    EXPECT_TRUE((std::ranges::ForwardRange<std::forward_list<char>>)); // `.size()` missing
    EXPECT_TRUE((std::ranges::ForwardRange<std::vector<char>>));
    EXPECT_TRUE((std::ranges::ForwardRange<std::deque<char>>));
    EXPECT_TRUE((std::ranges::ForwardRange<std::string>));

    EXPECT_TRUE((std::ranges::ForwardRange<concatenated_sequences_string_t>));
    EXPECT_TRUE((std::ranges::ForwardRange<seqan3::concatenated_sequences<std::vector<char>>>));
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
    // SequenceContainer_modified_by_const_iterator_bug is no issue anymore
    // it =
    s.insert(s.begin() + s.find_first_of(':') + 1, 2, '=');
    EXPECT_EQ("Exemplar is an:== example", s);

    // insert(const_iterator pos, InputIt first, InputIt last)
    {
        string_t seq = " string";
        //TODO should return type::iterator on strings, remove if
        // SequenceContainer_modified_by_const_iterator_bug is no issue anymore
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
    // concept bool SequenceContainer = requires (type val, type val2)
    //              ^~~~~~~~~~~~~~~~
    // /include/seqan3/range/container/concept.hpp:113:14: note:     with ‘std::basic_string<char> val’
    // /include/seqan3/range/container/concept.hpp:113:14: note:     with ‘std::basic_string<char> val2’
    // [...]
    // /include/seqan3/range/container/concept.hpp:113:14: note: the required expression ‘val.erase(val.cbegin(), val.cend())’ would be ill-formed

    using namespace std::string_literals;

    static_assert(detail::SequenceContainer_modified_by_const_iterator<string_t>);
    static_assert(!detail::SequenceContainer_modified_by_const_iterator_bug<string_t> && std::is_same_v<string_t, std::string>);
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

    if constexpr(!detail::SequenceContainer_modified_by_const_iterator_bug<>)
        container_concept_const_travis_bug_test<std::string>();
}

TEST(container_concept, SequenceContainer)
{
    EXPECT_FALSE((seqan3::SequenceContainer<std::array<char, 2>>));
    EXPECT_TRUE((seqan3::SequenceContainer<std::list<char>>));
    EXPECT_FALSE((seqan3::SequenceContainer<std::forward_list<char>>));
    EXPECT_TRUE((seqan3::SequenceContainer<std::vector<char>>));
    EXPECT_TRUE((seqan3::SequenceContainer<std::deque<char>>));
    EXPECT_TRUE((seqan3::SequenceContainer<std::string>));

    EXPECT_TRUE((seqan3::SequenceContainer<concatenated_sequences_string_t>));
    EXPECT_TRUE((seqan3::SequenceContainer<seqan3::concatenated_sequences<std::vector<char>>>));
}

TEST(container_concept, SequenceContainer_modified_by_const_iterator)
{
    EXPECT_FALSE((seqan3::detail::SequenceContainer_modified_by_const_iterator<std::array<char, 2>>));
    EXPECT_TRUE((seqan3::detail::SequenceContainer_modified_by_const_iterator<std::list<char>>));
    EXPECT_FALSE((seqan3::detail::SequenceContainer_modified_by_const_iterator<std::forward_list<char>>));
    EXPECT_TRUE((seqan3::detail::SequenceContainer_modified_by_const_iterator<std::vector<char>>));
    EXPECT_TRUE((seqan3::detail::SequenceContainer_modified_by_const_iterator<std::deque<char>>));
    if constexpr(!detail::SequenceContainer_modified_by_const_iterator_bug<>)
    {
        EXPECT_TRUE((seqan3::detail::SequenceContainer_modified_by_const_iterator<std::string>));
    }

    EXPECT_TRUE((seqan3::detail::SequenceContainer_modified_by_const_iterator<concatenated_sequences_string_t>));
    EXPECT_TRUE((seqan3::detail::SequenceContainer_modified_by_const_iterator<seqan3::concatenated_sequences<std::vector<char>>>));
}

TEST(container_concept, RandomAccessContainer)
{
    EXPECT_FALSE((seqan3::RandomAccessContainer<std::array<char, 2>>));
    EXPECT_FALSE((seqan3::RandomAccessContainer<std::list<char>>));
    EXPECT_FALSE((seqan3::RandomAccessContainer<std::forward_list<char>>));
    EXPECT_TRUE((seqan3::RandomAccessContainer<std::vector<char>>));
    EXPECT_TRUE((seqan3::RandomAccessContainer<std::deque<char>>));
    EXPECT_TRUE((seqan3::RandomAccessContainer<std::string>));

    EXPECT_TRUE((seqan3::RandomAccessContainer<concatenated_sequences_string_t>));
    EXPECT_TRUE((seqan3::RandomAccessContainer<seqan3::concatenated_sequences<std::vector<char>>>));
}

TEST(container_concept, ReservableContainer)
{
    EXPECT_FALSE((seqan3::ReservableContainer<std::array<char, 2>>));
    EXPECT_FALSE((seqan3::ReservableContainer<std::list<char>>));
    EXPECT_FALSE((seqan3::ReservableContainer<std::forward_list<char>>));
    EXPECT_TRUE((seqan3::ReservableContainer<std::vector<char>>));
    EXPECT_FALSE((seqan3::ReservableContainer<std::deque<char>>));
    EXPECT_TRUE((seqan3::ReservableContainer<std::string>));

    EXPECT_TRUE((seqan3::ReservableContainer<concatenated_sequences_string_t>));
    EXPECT_TRUE((seqan3::ReservableContainer<seqan3::concatenated_sequences<std::vector<char>>>));

    EXPECT_TRUE((seqan3::ReservableContainer<sdsl::bit_vector>));
    EXPECT_TRUE((seqan3::ReservableContainer<sdsl::int_vector<>>));
    EXPECT_TRUE((seqan3::ReservableContainer<sdsl::int_vector<13>>));
    EXPECT_TRUE((seqan3::ReservableContainer<sdsl::int_vector<64>>));
}

/* Check the SDSL containers */
//TODO

/* Check our containers */
//TODO
