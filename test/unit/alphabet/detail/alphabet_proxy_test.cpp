// SPDX-FileCopyrightText: 2006-2025 Knut Reinert & Freie Universität Berlin
// SPDX-FileCopyrightText: 2016-2025 Knut Reinert & MPI für molekulare Genetik
// SPDX-License-Identifier: BSD-3-Clause

#include <seqan3/alphabet/detail/alphabet_proxy.hpp>
#include <seqan3/alphabet/nucleotide/dna4.hpp>

#include "../alphabet_constexpr_test_template.hpp"
#include "../alphabet_test_template.hpp"
#include "../semi_alphabet_constexpr_test_template.hpp"
#include "../semi_alphabet_test_template.hpp"

class alphabet_proxy_example : public seqan3::alphabet_proxy<alphabet_proxy_example, seqan3::dna4>
{
private:
    using alphabet_type = seqan3::dna4;
    using base_t = seqan3::alphabet_proxy<alphabet_proxy_example, alphabet_type>;
    friend base_t;

    constexpr void on_update() noexcept
    {}

public:
    constexpr alphabet_proxy_example() noexcept = default;
    constexpr alphabet_proxy_example(alphabet_proxy_example const &) = default;
    constexpr alphabet_proxy_example(alphabet_proxy_example &&) = default;
    constexpr alphabet_proxy_example & operator=(alphabet_proxy_example const &) = default;
    constexpr alphabet_proxy_example & operator=(alphabet_proxy_example &&) = default;
    ~alphabet_proxy_example() = default;

    constexpr alphabet_proxy_example(alphabet_type const a) noexcept : base_t{a} {};

    using base_t::operator=;
};

INSTANTIATE_TYPED_TEST_SUITE_P(alphabet_proxy, alphabet, alphabet_proxy_example, );
INSTANTIATE_TYPED_TEST_SUITE_P(alphabet_proxy, semi_alphabet_test, alphabet_proxy_example, );
INSTANTIATE_TYPED_TEST_SUITE_P(alphabet_proxy, alphabet_constexpr, alphabet_proxy_example, );
INSTANTIATE_TYPED_TEST_SUITE_P(alphabet_proxy, semi_alphabet_constexpr, alphabet_proxy_example, );

// -----------------------------------------------------------------------------------------------------
// check handling of external types that do not provide members
// -----------------------------------------------------------------------------------------------------

namespace my_namespace
{

class my_alph
{
public:
    bool rank{};

    constexpr my_alph() noexcept = default;
    constexpr my_alph(my_alph const &) = default;
    constexpr my_alph & operator=(my_alph const &) = default;

    constexpr my_alph(bool rank) : rank{rank}
    {}

    constexpr friend bool operator==(my_alph lhs, my_alph rhs)
    {
        return lhs.rank == rhs.rank;
    }
    constexpr friend bool operator!=(my_alph lhs, my_alph rhs)
    {
        return lhs.rank != rhs.rank;
    }
    constexpr friend bool operator<=(my_alph lhs, my_alph rhs)
    {
        return lhs.rank <= rhs.rank;
    }
    constexpr friend bool operator>=(my_alph lhs, my_alph rhs)
    {
        return lhs.rank >= rhs.rank;
    }
    constexpr friend bool operator<(my_alph lhs, my_alph rhs)
    {
        return lhs.rank < rhs.rank;
    }
    constexpr friend bool operator>(my_alph lhs, my_alph rhs)
    {
        return lhs.rank > rhs.rank;
    }
};

constexpr size_t alphabet_size(my_alph const &) noexcept
{
    return 2;
}

constexpr bool to_rank(my_alph const a) noexcept
{
    return a.rank;
}

constexpr my_alph & assign_rank_to(bool const r, my_alph & a) noexcept
{
    a.rank = r;
    return a;
}

constexpr char to_char(my_alph const a) noexcept
{
    if (a.rank)
        return '1';
    else
        return '0';
}

constexpr my_alph & assign_char_to(char const c, my_alph & a) noexcept
{
    switch (c)
    {
    case '0':
    case 'F':
    case 'f':
        a.rank = 0;
        return a;
    default:
        a.rank = 1;
        return a;
    }
}

} // namespace my_namespace

static_assert(seqan3::alphabet_size<my_namespace::my_alph> == 2);
static_assert(seqan3::semialphabet<my_namespace::my_alph>);
static_assert(seqan3::alphabet<my_namespace::my_alph>);

class alphabet_proxy_example2 : public seqan3::alphabet_proxy<alphabet_proxy_example2, my_namespace::my_alph>
{
private:
    using alphabet_type = my_namespace::my_alph;
    using base_t = seqan3::alphabet_proxy<alphabet_proxy_example2, alphabet_type>;
    friend base_t;

    constexpr void on_update() noexcept
    {}

public:
    constexpr alphabet_proxy_example2() noexcept = default;
    constexpr alphabet_proxy_example2(alphabet_proxy_example2 const &) = default;
    constexpr alphabet_proxy_example2(alphabet_proxy_example2 &&) = default;
    constexpr alphabet_proxy_example2 & operator=(alphabet_proxy_example2 const &) = default;
    constexpr alphabet_proxy_example2 & operator=(alphabet_proxy_example2 &&) = default;
    ~alphabet_proxy_example2() = default;

    constexpr alphabet_proxy_example2(alphabet_type const a) noexcept : base_t{a} {};
};

INSTANTIATE_TYPED_TEST_SUITE_P(alphabet_proxy2, alphabet, alphabet_proxy_example2, );
INSTANTIATE_TYPED_TEST_SUITE_P(alphabet_proxy2, alphabet_constexpr, alphabet_proxy_example2, );
