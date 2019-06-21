// -----------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2019, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2019, Knut Reinert & MPI für molekulare Genetik
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at: https://github.com/seqan/seqan3/blob/master/LICENSE
// -----------------------------------------------------------------------------------------------------

#include <seqan3/alphabet/concept.hpp>   // for seqan3::Alphabet concept checks

struct dna2
{
    uint8_t rank{};

//! [writable_semialphabet]
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
//! [writable_semialphabet]

//! [alphabet]
    char to_char() const noexcept
    {
        char const rank_to_char[2] {'S', 'W'};
        return rank_to_char[rank];
    }
//! [alphabet]

//! [writable_alphabet]
    dna2 & assign_char(char const ch) noexcept
    {
        switch (ch)
        {
            case 'W': rank = 1; break;
            default:  rank = 0;
        }
        return *this;
    }

    static bool char_is_valid(char const ch) noexcept
    {
        return (ch == dna2{}.assign_char(ch).to_char());
    }
//! [writable_alphabet]
};

//! [equality]
bool operator==(dna2 const & lhs, dna2 const & rhs)
{
    return lhs.rank == rhs.rank;
}
//! [equality]

//! [inequality]
bool operator!=(dna2 const & lhs, dna2 const & rhs)
{
    return !(lhs == rhs);
}
//! [inequality]

//! [comparison]
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

static_assert(std::StrictTotallyOrdered<dna2>);   // ok
//! [comparison]

//! [equality_comparable]
static_assert(std::EqualityComparable<dna2>);     // ok
static_assert(std::StrictTotallyOrdered<dna2>);   // compiler error
//! [equality_comparable]

static_assert(seqan3::Semialphabet<dna2>);
static_assert(std::CopyConstructible<dna2>);

//! [writable_alphabet_concept]
static_assert(seqan3::Alphabet<dna2>);            // ok
static_assert(seqan3::WritableAlphabet<dna2>);    // ok
//! [writable_alphabet_concept]
