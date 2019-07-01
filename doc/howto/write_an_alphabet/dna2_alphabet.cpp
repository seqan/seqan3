// -----------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2019, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2019, Knut Reinert & MPI für molekulare Genetik
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at: https://github.com/seqan/seqan3/blob/master/LICENSE
// -----------------------------------------------------------------------------------------------------

//! [writable_alphabet]
#include <seqan3/alphabet/concept.hpp>                   // alphabet concept checks

struct dna2
{
    uint8_t rank{};

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

    char to_char() const noexcept
    {
        // map 0 => 'S' and 1 => 'W'
        char const rank_to_char[2] {'S', 'W'};
        return rank_to_char[rank];
    }

    dna2 & assign_char(char const ch) noexcept
    {
        switch (ch)
        {
            case 'W': rank = 1;
            case 'w': rank = 1; break;                   // allow assignment from uppercase and lowercase
            default:  rank = 0;                          // unknown characters are mapped to 0 (=> 'S')
        }
        return *this;
    }

    static bool char_is_valid(char const ch) noexcept
    {
        return (ch == dna2{}.assign_char(ch).to_char());
    }
};

bool operator==(dna2 const & lhs, dna2 const & rhs)
{
    return lhs.rank == rhs.rank;
}

bool operator!=(dna2 const & lhs, dna2 const & rhs)
{
    return !(lhs == rhs);
}

bool operator<(dna2 const & lhs, dna2 const & rhs)
{
    return lhs.rank < rhs.rank;
}

bool operator<=(dna2 const & lhs, dna2 const & rhs)
{
    return lhs.rank <= rhs.rank;
}

bool operator>(dna2 const & lhs, dna2 const & rhs)
{
    return lhs.rank > rhs.rank;
}

bool operator>=(dna2 const & lhs, dna2 const & rhs)
{
    return lhs.rank >= rhs.rank;
}
//! [writable_alphabet]

//! [writable_alphabet_concept]
static_assert(seqan3::Alphabet<dna2>);                   // ok
static_assert(seqan3::WritableAlphabet<dna2>);           // ok
//! [writable_alphabet_concept]

//! [dummy_requirement]
template <seqan3::Alphabet check_this_type>
void foo()
{}

int main()
{
    foo<dna2>();
}
//! [dummy_requirement]
