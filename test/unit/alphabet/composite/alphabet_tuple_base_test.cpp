// SPDX-FileCopyrightText: 2006-2025 Knut Reinert & Freie Universität Berlin
// SPDX-FileCopyrightText: 2016-2025 Knut Reinert & MPI für molekulare Genetik
// SPDX-License-Identifier: BSD-3-Clause

#include <gtest/gtest.h>

#include <seqan3/alphabet/composite/alphabet_tuple_base.hpp>
#include <seqan3/alphabet/nucleotide/rna4.hpp>
#include <seqan3/alphabet/nucleotide/rna5.hpp>
#include <seqan3/core/debug_stream/tuple.hpp>

#include "../semi_alphabet_test_template.hpp"
#include "alphabet_tuple_base_test_template.hpp"

using seqan3::operator""_dna4;
using seqan3::operator""_dna5;
using seqan3::operator""_rna4;
using seqan3::operator""_rna5;

template <typename type1, typename type2>
struct test_composite : public seqan3::alphabet_tuple_base<test_composite<type1, type2>, type1, type2>
{
    using base_t = seqan3::alphabet_tuple_base<test_composite<type1, type2>, type1, type2>;
    using base_t::base_t;
    using base_t::operator=;
};

template <>
class alphabet_tuple_base_test<test_composite<seqan3::dna4, seqan3::dna5>> : public ::testing::Test
{
public:
    using T = test_composite<seqan3::dna4, seqan3::dna5>;

    T instance = T{value_1(), value_2()};
    T zero_instance = T{decltype(value_1()){}, decltype(value_2()){}};
    size_t tup_size{2};

    // test_composite<seqan3::dna4, seqan3::dna5>
    // -------------------------------------------------------------------------
    seqan3::dna4 value_1()
    {
        return 'G'_dna4;
    }
    seqan3::rna4 assignable_to_value_1()
    {
        return 'G'_rna4;
    }
    seqan3::dna5 value_2()
    {
        return 'G'_dna5;
    }
    seqan3::rna5 assignable_to_value_2()
    {
        return 'G'_rna5;
    }
    auto values_to_cmp()
    {
        return std::make_tuple(/*low */ 'A'_dna4,
                               'A'_dna5,
                               /*mid */ 'C'_dna4,
                               'C'_dna5,
                               /*high*/ 'T'_dna4,
                               'T'_dna5);
    }
};

using test_composite_types = ::testing::Types<test_composite<seqan3::dna4, seqan3::dna5>>;

INSTANTIATE_TYPED_TEST_SUITE_P(test_composite, semi_alphabet_test, test_composite_types, );
INSTANTIATE_TYPED_TEST_SUITE_P(test_composite, alphabet_tuple_base_test, test_composite_types, );
