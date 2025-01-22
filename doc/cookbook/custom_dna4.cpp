// SPDX-FileCopyrightText: 2006-2025 Knut Reinert & Freie Universität Berlin
// SPDX-FileCopyrightText: 2016-2025 Knut Reinert & MPI für molekulare Genetik
// SPDX-License-Identifier: CC0-1.0

#include <seqan3/alphabet/nucleotide/nucleotide_base.hpp>
#include <seqan3/core/debug_stream.hpp>

// We inherit from seqan3::nucleotide_base s.t. we do not need to implement the full nucleotide interface
// but it is sufficient to define `rank_to_char`, `char_to_rank`, and `complement_table`.
class my_dna4 : public seqan3::nucleotide_base<my_dna4, 4 /*alphabet size is 4*/>
{
public:
    using nucleotide_base<my_dna4, 4>::nucleotide_base; // Use constructors of the base class.

private:
    // Returns the character representation of rank. This is where rank conversion for to_char() is handled!
    static constexpr char_type rank_to_char(rank_type const rank)
    {
        return rank_to_char_table[rank];
    }

    // Returns the rank representation of character. This is where char conversion for assign_char() is handled!
    static constexpr rank_type char_to_rank(char_type const chr)
    {
        using index_t = std::make_unsigned_t<char_type>;
        return char_to_rank_table[static_cast<index_t>(chr)];
    }

    // Returns the complement by rank. This is where complement is handled and with this, my_dna4 models
    // seqan3::nucleotide_alphabet.
    static constexpr rank_type rank_complement(rank_type const rank)
    {
        return rank_complement_table[rank];
    }

private:
    // === lookup-table implementation detail ===

    // Value to char conversion table.
    static constexpr char_type rank_to_char_table[alphabet_size]{'A', 'C', 'G', 'T'}; // rank 0,1,2,3

    // Char-to-value conversion table.
    static constexpr std::array<rank_type, 256> char_to_rank_table{
        []() constexpr
        {
            // By default, everything has rank 0 which equals `A`.
            std::array<rank_type, 256> conversion_table{};

            conversion_table['C'] = conversion_table['c'] = 1;
            conversion_table['G'] = conversion_table['g'] = 2;
            conversion_table['T'] = conversion_table['t'] = 3;
            conversion_table['U'] = conversion_table['T']; // set U equal to T
            conversion_table['u'] = conversion_table['t']; // set u equal to t

            return conversion_table;
        }()};

    // The rank complement table.
    static constexpr rank_type rank_complement_table[alphabet_size]{
        3, // T is complement of 'A'_dna4
        2, // G is complement of 'C'_dna4
        1, // C is complement of 'G'_dna4
        0  // A is complement of 'T'_dna4
    };

    friend nucleotide_base<my_dna4, 4>;         // Grant seqan3::nucleotide_base access to private/protected members.
    friend nucleotide_base<my_dna4, 4>::base_t; // Grant seqan3::alphabet_base access to private/protected members.
};

// Defines the `_my_dna4` *char literal* so you can write `'C'_my_dna4` instead of `my_dna4{}.assign_char('C')`.
constexpr my_dna4 operator""_my_dna4(char const c) noexcept
{
    return my_dna4{}.assign_char(c);
}

int main()
{
    my_dna4 my_letter{'C'_my_dna4};

    my_letter.assign_char('S'); // Characters other than A,C,G,T are implicitly converted to `A`.

    seqan3::debug_stream << my_letter << "\n";                     // "A";
    seqan3::debug_stream << seqan3::complement(my_letter) << "\n"; // "T";
}
