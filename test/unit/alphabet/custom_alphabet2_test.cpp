// SPDX-FileCopyrightText: 2006-2025 Knut Reinert & Freie Universität Berlin
// SPDX-FileCopyrightText: 2016-2025 Knut Reinert & MPI für molekulare Genetik
// SPDX-License-Identifier: BSD-3-Clause

#include <seqan3/core/platform.hpp>

// Tests the ADL-capabilities of the alphabet customisation point objects; type that isn't default constructible

//![my_alph]
#include <cstddef>     // for size_t
#include <type_traits> // for std::type_identity

#include <seqan3/alphabet/concept.hpp> // for seqan3::alphabet

namespace my_namespace
{

class my_alph
{
public:
    bool rank;

    my_alph() = delete;
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

constexpr size_t alphabet_size(std::type_identity<my_alph> const &) noexcept
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

constexpr bool char_is_valid_for(char const c, std::type_identity<my_alph> const &) noexcept
{
    switch (c)
    {
    case '0':
    case 'F':
    case 'f':
    case '1':
    case 'T':
    case 't':
        return true;
    default:
        return false;
    }
}

} // namespace my_namespace

static_assert(seqan3::alphabet_size<my_namespace::my_alph> == 2);
static_assert(seqan3::char_is_valid_for<my_namespace::my_alph>('T'));
static_assert(!seqan3::char_is_valid_for<my_namespace::my_alph>('!'));
static_assert(seqan3::semialphabet<my_namespace::my_alph>);
static_assert(seqan3::alphabet<my_namespace::my_alph>);
//![my_alph]

// Not tested with rest of test-suite because the test-suite relies on default-constructibility
