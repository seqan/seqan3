// -----------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2019, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2019, Knut Reinert & MPI für molekulare Genetik
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at: https://github.com/seqan/seqan3/blob/master/LICENSE
// -----------------------------------------------------------------------------------------------------

//! [exercise]
#include <iostream>   // std::cerr, std::endl
#include <seqan3/alphabet/all.hpp>

using namespace seqan3;

class dna2
{
public:
    using rank_type = uint8_t;
    using char_type = char;

    static rank_type const value_size = 2;

    rank_type to_rank() const noexcept
    {
        return rank;
    };

    char_type to_char() const noexcept
    {
        char_type const rank_to_char[2] {'S', 'W'};
        return rank_to_char[rank];
    }

    dna2 & assign_rank(rank_type const rk) noexcept
    {
        assert(rk < value_size);
        rank = rk;
        return *this;
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
            throw invalid_char_assignment{"dna2", ch};

        return assign_char(ch);
    }

private:
    rank_type rank{};
};

bool operator==(dna2 const & lhs, dna2 const & rhs)
{
    return lhs.to_rank() == rhs.to_rank();
}

bool operator!=(dna2 const & lhs, dna2 const & rhs)
{
    return lhs.to_rank() != rhs.to_rank();
}

bool operator<(dna2 const & lhs, dna2 const & rhs)
{
    return lhs.to_rank() < rhs.to_rank();
}

bool operator>(dna2 const & lhs, dna2 const & rhs)
{
    return lhs.to_rank() > rhs.to_rank();
}

bool operator<=(dna2 const & lhs, dna2 const & rhs)
{
    return lhs.to_rank() <= rhs.to_rank();
}

bool operator>=(dna2 const & lhs, dna2 const & rhs)
{
    return lhs.to_rank() >= rhs.to_rank();
}

// Constrained function that works only for seqan3::Alphabet types.
template <typename alph_type>
    requires Alphabet<alph_type>
void test_function(alph_type)
{
    std::cerr << "You're good!" << std::endl;
    std::cerr << "The alphabet size is " << (unsigned)alphabet_size_v<dna2> << "." << std::endl;
}

int main ()
{
    // Let's test our new alphabet class here. The compilation fails, if members are missing.
    test_function(dna2{});
    return 0;
}
//! [exercise]
