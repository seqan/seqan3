// SPDX-FileCopyrightText: 2006-2025 Knut Reinert & Freie Universität Berlin
// SPDX-FileCopyrightText: 2016-2025 Knut Reinert & MPI für molekulare Genetik
// SPDX-License-Identifier: BSD-3-Clause

#include <seqan3/core/platform.hpp>

#include "alphabet_constexpr_test_template.hpp"
#include "alphabet_test_template.hpp"
#include "semi_alphabet_constexpr_test_template.hpp"
#include "semi_alphabet_test_template.hpp"

// Tests the capabilities of the explicit alphabet customisation

//![third_party_type]
#include <cstddef> // for size_t

#include <seqan3/alphabet/concept.hpp> // for seqan3::alphabet

// this is from some other library:
namespace third_party_ns
{

enum class third_party_type
{
    zero,
    one,
    two
};

} // namespace third_party_ns

// ------------------------------------------------------------------------------------

// this is in your code (no namespace):
template <>
struct seqan3::custom::alphabet<third_party_ns::third_party_type>
{
    using alphabet_t = third_party_ns::third_party_type;

    static constexpr size_t alphabet_size = 3;

    static constexpr size_t to_rank(alphabet_t const a) noexcept
    {
        return static_cast<size_t>(a);
    }

    static constexpr alphabet_t & assign_rank_to(size_t const r, alphabet_t & a) noexcept
    {
        switch (r)
        {
        case 0:
            a = alphabet_t::zero;
            return a;
        case 1:
            a = alphabet_t::one;
            return a;
        default:
            a = alphabet_t::two;
            return a;
        }
    }

    static constexpr char to_char(alphabet_t const a) noexcept
    {
        switch (a)
        {
        case alphabet_t::zero:
            return '0';
        case alphabet_t::one:
            return '1';
        default:
            return '2';
        }
    }

    static constexpr alphabet_t & assign_char_to(char const c, alphabet_t & a) noexcept
    {
        switch (c)
        {
        case '0':
            a = alphabet_t::zero;
            return a;
        case '1':
            a = alphabet_t::one;
            return a;
        default:
            a = alphabet_t::two;
            return a;
        }
    }
};

static_assert(seqan3::alphabet<third_party_ns::third_party_type>);
//![third_party_type]

INSTANTIATE_TYPED_TEST_SUITE_P(third_party_type, alphabet, third_party_ns::third_party_type, );
INSTANTIATE_TYPED_TEST_SUITE_P(third_party_type, semi_alphabet_test, third_party_ns::third_party_type, );
INSTANTIATE_TYPED_TEST_SUITE_P(third_party_type, alphabet_constexpr, third_party_ns::third_party_type, );
INSTANTIATE_TYPED_TEST_SUITE_P(third_party_type, semi_alphabet_constexpr, third_party_ns::third_party_type, );
