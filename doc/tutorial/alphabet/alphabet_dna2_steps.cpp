// -----------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2019, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2019, Knut Reinert & MPI für molekulare Genetik
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at: https://github.com/seqan/seqan3/blob/master/LICENSE
// -----------------------------------------------------------------------------------------------------

#include <seqan3/alphabet/concept.hpp>   // for seqan3::Alphabet concept checks
#include <seqan3/alphabet/exception.hpp> // for seqan3::invalid_char_assignment

struct dna2
{
    using rank_type = uint8_t;
    rank_type rank{};

//! [semialphabet]
    static rank_type const value_size = 2;

    rank_type to_rank() const noexcept
    {
        return rank;
    }

    dna2 & assign_rank(rank_type const rk) noexcept
    {
        assert(rk < value_size);
        rank = rk;
        return *this;
    }
//! [semialphabet]

//! [alphabet]
    using char_type = char;

    char_type to_char() const noexcept
    {
        char_type const rank_to_char[2] {'S', 'W'};
        return rank_to_char[rank];
    }

    dna2 & assign_char(char_type const ch) noexcept
    {
        switch (ch)
        {
            case 'W': rank = 1; break;
            default:  rank = 0;
        }
        return *this;
    }

    static bool char_is_valid(char_type const ch) noexcept
    {
        return (ch == dna2{}.assign_char(ch).to_char());
    }

    dna2 & assign_char_strict(char_type const ch)
    {
        if (!char_is_valid(ch))
            throw seqan3::invalid_char_assignment{"dna2", ch};

        return assign_char(ch);
    }
//! [alphabet]
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

//! [regular]
static_assert(std::Regular<dna2>); // ok
//! [regular]

//! [compare]
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

static_assert(std::StrictTotallyOrdered<dna2>); // ok
//! [compare]

static_assert(seqan3::Semialphabet<dna2>);

//! [alphabet_concept]
static_assert(seqan3::Alphabet<dna2>); // ok
//! [alphabet_concept]
