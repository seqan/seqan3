// -----------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2020, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2020, Knut Reinert & MPI für molekulare Genetik
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at: https://github.com/seqan/seqan3/blob/master/LICENSE.md
// -----------------------------------------------------------------------------------------------------

#include <seqan3/core/platform.hpp>

#include "alphabet_test_template.hpp"
#include "alphabet_constexpr_test_template.hpp"
#include "semi_alphabet_test_template.hpp"
#include "semi_alphabet_constexpr_test_template.hpp"

// Tests the ADL-capabilities of the alphabet customisation point objects

//![my_alph]
#include <cstddef>                      // for size_t
#include <seqan3/alphabet/concept.hpp>  // for seqan3::alphabet

namespace my_namespace
{

enum class my_alph
{
    ZERO,
    ONE,
    TWO
};

constexpr size_t alphabet_size(my_alph const &) noexcept
{
    return 3;
}

constexpr size_t to_rank(my_alph const a) noexcept
{
    return static_cast<size_t>(a);
}

constexpr my_alph & assign_rank_to(size_t const r, my_alph & a) noexcept
{
    switch (r)
    {
        case 0:  a = my_alph::ZERO; return a;
        case 1:  a = my_alph::ONE;  return a;
        default: a = my_alph::TWO;  return a;
    }
}

constexpr char to_char(my_alph const a) noexcept
{
    switch (a)
    {
        case my_alph::ZERO: return '0';
        case my_alph::ONE:  return '1';
        default:            return '2';
    }
}

constexpr my_alph & assign_char_to(char const c, my_alph & a) noexcept
{
    switch (c)
    {
        case '0': a = my_alph::ZERO; return a;
        case '1': a = my_alph::ONE;  return a;
        default:  a = my_alph::TWO;  return a;
    }
}

} // namespace my_namespace

static_assert(seqan3::alphabet<my_namespace::my_alph>);
//![my_alph]

INSTANTIATE_TYPED_TEST_SUITE_P(my_alph, alphabet_, my_namespace::my_alph, );
INSTANTIATE_TYPED_TEST_SUITE_P(my_alph, semi_alphabet_test, my_namespace::my_alph, );
INSTANTIATE_TYPED_TEST_SUITE_P(my_alph, alphabet_constexpr, my_namespace::my_alph, );
INSTANTIATE_TYPED_TEST_SUITE_P(my_alph, semi_alphabet_constexpr, my_namespace::my_alph, );
