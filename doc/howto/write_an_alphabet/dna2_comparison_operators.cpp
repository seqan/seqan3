// -----------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2019, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2019, Knut Reinert & MPI für molekulare Genetik
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at: https://github.com/seqan/seqan3/blob/master/LICENSE
// -----------------------------------------------------------------------------------------------------

//! [comparison]
#include <seqan3/alphabet/concept.hpp>                   // alphabet concept checks

struct dna2
{
    uint8_t rank{};
};

// Equality and inequality operators

bool operator==(dna2 const & lhs, dna2 const & rhs)
{
    return lhs.rank == rhs.rank;
}

bool operator!=(dna2 const & lhs, dna2 const & rhs)
{
    return !(lhs == rhs);
}

// Comparison operators

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

static_assert(std::StrictTotallyOrdered<dna2>);          // ok
//! [comparison]

//! [equality_comparable_concept]
static_assert(std::EqualityComparable<dna2>);            // ok
//! [equality_comparable_concept]
