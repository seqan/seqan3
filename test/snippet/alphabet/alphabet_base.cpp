// SPDX-FileCopyrightText: 2006-2025 Knut Reinert & Freie Universität Berlin
// SPDX-FileCopyrightText: 2016-2025 Knut Reinert & MPI für molekulare Genetik
// SPDX-License-Identifier: CC0-1.0

#include <seqan3/alphabet/alphabet_base.hpp>
#include <seqan3/alphabet/concept.hpp>
#include <seqan3/utility/char_operations/transform.hpp> // seqan3::to_lower

class ab : public seqan3::alphabet_base<ab, 2>
{
private:
    // make the base class a friend so it can access the tables:
    friend alphabet_base<ab, 2>;

    // This function is expected by seqan3::alphabet_base
    static constexpr char_type rank_to_char(rank_type const rank)
    {
        // via a lookup table
        return rank_to_char_table[rank];
        // or via an arithmetic expression
        return rank == 1 ? 'B' : 'A';
    }

    // This function is expected by seqan3::alphabet_base
    static constexpr rank_type char_to_rank(char_type const chr)
    {
        // via a lookup table
        using index_t = std::make_unsigned_t<char_type>;
        return char_to_rank_table[static_cast<index_t>(chr)];
        // or via an arithmetic expression
        return seqan3::to_lower(chr) == 'b' ? 1 : 0;
    }

private:
    // === lookup-table implementation detail ===

    // map 0 -> A and 1 -> B
    static constexpr std::array<char_type, alphabet_size> rank_to_char_table{'A', 'B'};

    // map every letter to rank zero, except Bs
    static constexpr std::array<rank_type, 256> char_to_rank_table{
        // initialise with an immediately evaluated lambda expression:
        []()
        {
            std::array<rank_type, 256> ret{}; // initialise all values with 0 / 'A'

            // only 'b' and 'B' result in rank 1
            ret['b'] = 1;
            ret['B'] = 1;

            return ret;
        }()};
};

// The class ab satisfies the alphabet concept.
static_assert(seqan3::alphabet<ab>);
static_assert(seqan3::writable_alphabet<ab>);
