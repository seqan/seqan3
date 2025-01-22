// SPDX-FileCopyrightText: 2006-2025 Knut Reinert & Freie Universität Berlin
// SPDX-FileCopyrightText: 2016-2025 Knut Reinert & MPI für molekulare Genetik
// SPDX-License-Identifier: BSD-3-Clause

//! [writable_alphabet]
#include <cassert>

#include <seqan3/alphabet/concept.hpp> // alphabet concept checks

struct dna2
{
    uint8_t rank{};

    // semialphabet

    static constexpr size_t alphabet_size = 2;

    uint8_t to_rank() const noexcept
    {
        return rank;
    }

    dna2 & assign_rank(uint8_t const rk) noexcept
    {
        assert(rk < alphabet_size);
        rank = rk;
        return *this;
    }

    // alphabet

    char to_char() const noexcept
    {
        // map 0 => 'S' and 1 => 'W'
        char const rank_to_char[2]{'S', 'W'};
        return rank_to_char[rank];
    }

    dna2 & assign_char(char const ch) noexcept
    {
        switch (ch)
        {
        case 'W':
        case 'w':
            rank = 1;
            break; // allow assignment from uppercase and lowercase
        default:
            rank = 0; // unknown characters are mapped to 0 (=> 'S')
        }
        return *this;
    }

    // Optional: Can be omitted.

    static bool char_is_valid(char const ch) noexcept
    {
        return (ch == dna2{}.assign_char(ch).to_char());
    }

    // Equality and inequality operators

    friend bool operator==(dna2 const & lhs, dna2 const & rhs) noexcept
    {
        return lhs.rank == rhs.rank;
    }

    friend bool operator!=(dna2 const & lhs, dna2 const & rhs) noexcept
    {
        return !(lhs == rhs);
    }

    // Comparison operators

    friend bool operator<(dna2 const & lhs, dna2 const & rhs) noexcept
    {
        return lhs.rank < rhs.rank;
    }

    friend bool operator<=(dna2 const & lhs, dna2 const & rhs) noexcept
    {
        return lhs.rank <= rhs.rank;
    }

    friend bool operator>(dna2 const & lhs, dna2 const & rhs) noexcept
    {
        return lhs.rank > rhs.rank;
    }

    friend bool operator>=(dna2 const & lhs, dna2 const & rhs) noexcept
    {
        return lhs.rank >= rhs.rank;
    }
};
//! [writable_alphabet]

//! [writable_alphabet_concept]
static_assert(seqan3::alphabet<dna2>);          // ok
static_assert(seqan3::writable_alphabet<dna2>); // ok
//! [writable_alphabet_concept]

//! [dummy_requirement]
template <seqan3::alphabet check_this_type>
void foo()
{}

int main()
{
    foo<dna2>();
}
//! [dummy_requirement]
