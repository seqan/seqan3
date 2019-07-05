// -----------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2019, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2019, Knut Reinert & MPI für molekulare Genetik
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at: https://github.com/seqan/seqan3/blob/master/LICENSE
// -----------------------------------------------------------------------------------------------------

//! [semialphabet]
#include <seqan3/alphabet/concept.hpp>                   // alphabet concept checks

struct dna2
{
    uint8_t rank{};

    // Semialphabet

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

    // Equality and inequality operators ...
    // Comparison operators ...
};
//! [semialphabet]

bool operator==(dna2 const & lhs, dna2 const & rhs) noexcept
{
    return lhs.rank == rhs.rank;
}

bool operator!=(dna2 const & lhs, dna2 const & rhs) noexcept
{
    return !(lhs == rhs);
}

bool operator<(dna2 const & lhs, dna2 const & rhs) noexcept
{
    return lhs.rank < rhs.rank;
}

bool operator<=(dna2 const & lhs, dna2 const & rhs) noexcept
{
    return lhs.rank <= rhs.rank;
}

bool operator>(dna2 const & lhs, dna2 const & rhs) noexcept
{
    return lhs.rank > rhs.rank;
}

bool operator>=(dna2 const & lhs, dna2 const & rhs) noexcept
{
    return lhs.rank >= rhs.rank;
}

//! [writable_semialphabet_concept]
static_assert(seqan3::Semialphabet<dna2>);               // ok
static_assert(seqan3::WritableSemialphabet<dna2>);       // ok
//! [writable_semialphabet_concept]

//! [free_functions]
int main ()
{
    dna2 chr{};
    seqan3::assign_rank_to(1, chr);                      // chr is assigned rank 1
    uint8_t rnk = seqan3::to_rank(chr);                  // query rank value  => 1
}
//! [free_functions]
