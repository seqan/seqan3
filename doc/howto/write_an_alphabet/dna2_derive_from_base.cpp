// -----------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2020, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2020, Knut Reinert & MPI für molekulare Genetik
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at: https://github.com/seqan/seqan3/blob/master/LICENSE
// -----------------------------------------------------------------------------------------------------

//! [dna2]
#include <array>                                         // std::array
#include <seqan3/alphabet/alphabet_base.hpp>             // alphabet_base
#include <seqan3/alphabet/concept.hpp>                   // alphabet concept checks

// derive from alphabet_base
struct dna2 : public seqan3::alphabet_base<dna2, 2>
{
private:
    // map 0 => 'S' and 1 => 'W'
    static constexpr char_type rank_to_char[alphabet_size] {'S', 'W'};

    static constexpr std::array<rank_type, 256> char_to_rank
    {
        // initialise with an immediately evaluated lambda expression:
        [] () constexpr
        {
            std::array<rank_type, 256> ret{}; // initialise all values with 0 (=> 'S')
            ret['W'] = 1;                     // only 'W' and 'w' result in rank 1
            ret['w'] = 1;
            return ret;
        } ()
    };

    // make the base class a friend so it can access the tables:
    friend seqan3::alphabet_base<dna2, 2>;
};

// check the concepts
static_assert(seqan3::alphabet<dna2>);                   // ok
static_assert(seqan3::writable_alphabet<dna2>);           // ok
//! [dna2]
