// -----------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2019, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2019, Knut Reinert & MPI für molekulare Genetik
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at: https://github.com/seqan/seqan3/blob/master/LICENSE.md
// -----------------------------------------------------------------------------------------------------

#include <gtest/gtest.h>

#include <seqan3/alphabet/composite/alphabet_tuple_base.hpp>
#include <seqan3/alphabet/nucleotide/rna4.hpp>
#include <seqan3/alphabet/nucleotide/rna5.hpp>

#include "alphabet_tuple_base_test_template.hpp"

using namespace seqan3;

template <typename type1, typename type2>
struct test_composite : public alphabet_tuple_base<test_composite<type1, type2>, type1, type2>
{
    using base_t = alphabet_tuple_base<test_composite<type1, type2>, type1, type2>;
    using base_t::base_t;
    using base_t::operator=;

    using base_t::operator==;
    using base_t::operator!=;
    using base_t::operator<;
    using base_t::operator>;
    using base_t::operator<=;
    using base_t::operator>=;
};

template <>
class alphabet_tuple_base_test<test_composite<dna4, dna5>> : public ::testing::Test
{
public:
    using T = test_composite<dna4, dna5>;

    T instance = T{value_1(), value_2()};
    T zero_instance = T{decltype(value_1()){}, decltype(value_2()){}};
    size_t tup_size{2};

    // test_composite<dna4, dna5>
    // -------------------------------------------------------------------------
    dna4 value_1()
    {
        return 'G'_dna4;
    }
    rna4 assignable_to_value_1()
    {
        return 'G'_rna4;
    }
    dna5 value_2()
    {
        return 'G'_dna5;
    }
    rna5 assignable_to_value_2()
    {
        return 'G'_rna5;
    }
    auto values_to_cmp()
    {
        return std::make_tuple(/*low */'A'_dna4, 'A'_dna5,
                               /*mid */'C'_dna4, 'C'_dna5,
                               /*high*/'T'_dna4, 'T'_dna5);
    }
};

using test_composite_types = ::testing::Types<test_composite<dna4, dna5>>;

INSTANTIATE_TYPED_TEST_CASE_P(test_composite, alphabet_tuple_base_test, test_composite_types);
