// SPDX-FileCopyrightText: 2006-2025 Knut Reinert & Freie Universität Berlin
// SPDX-FileCopyrightText: 2016-2025 Knut Reinert & MPI für molekulare Genetik
// SPDX-License-Identifier: CC0-1.0

#include <seqan3/alphabet/aminoacid/concept.hpp>

namespace your_namespace
{

// your own aminoacid definition
struct your_aa : seqan3::aminoacid_empty_base
{
    //...
};

} // namespace your_namespace

static_assert(seqan3::enable_aminoacid<your_namespace::your_aa> == true);

/***** OR *****/

namespace your_namespace2
{

// your own aminoacid definition
struct your_aa
{
    //...
};

constexpr bool enable_aminoacid(your_aa) noexcept
{
    return true;
}

} // namespace your_namespace2

static_assert(seqan3::enable_aminoacid<your_namespace2::your_aa> == true);
