// SPDX-FileCopyrightText: 2006-2025 Knut Reinert & Freie Universität Berlin
// SPDX-FileCopyrightText: 2016-2025 Knut Reinert & MPI für molekulare Genetik
// SPDX-License-Identifier: BSD-3-Clause

//! [equality]
#include <seqan3/alphabet/concept.hpp> // alphabet concept checks

struct dna2
{
    uint8_t rank{};

    // Equality operator

    friend bool operator==(dna2 const & lhs, dna2 const & rhs) noexcept
    {
        return lhs.rank == rhs.rank;
    }
};
//! [equality]
