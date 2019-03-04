// -----------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2019, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2019, Knut Reinert & MPI für molekulare Genetik
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at: https://github.com/seqan/seqan3/blob/master/LICENSE
// -----------------------------------------------------------------------------------------------------

//! [struct]
#include <seqan3/alphabet/concept.hpp>   // for seqan3::Alphabet concept checks
#include <seqan3/alphabet/exception.hpp> // for seqan3::invalid_char_assignment

struct dna2
{
    using rank_type = uint8_t;
    rank_type rank{};
};
//! [struct]
