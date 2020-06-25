#include <seqan3/alphabet/nucleotide/nucleotide_base.hpp>
#include <seqan3/core/debug_stream.hpp>

// We inherit from seqan3::nucleotide_base s.t. we do not need to implement the full nucleotide interface
// but it is sufficient to define `rank_to_char`, `char_to_rank`, and `complement_table`.
class my_dna4 : public seqan3::nucleotide_base<my_dna4, 4/*alphabet size is 4*/>
{
public:
    using nucleotide_base<my_dna4, 4>::nucleotide_base; // Use constructors of the base class.

private:
    // Value to char conversion table.
    static constexpr char_type rank_to_char[alphabet_size] {'A', 'C', 'G', 'T'}; // rank 0,1,2,3

    // Char-to-value conversion table. This is where char conversion for assign_char() is handled!
    static constexpr std::array<rank_type, 256> char_to_rank
    {
        [] () constexpr
        {
            std::array<rank_type, 256> conversion_table{}; // By default now, everything has rank 0 which equals `A`.

            conversion_table['C'] = conversion_table['c'] = 1;
            conversion_table['G'] = conversion_table['g'] = 2;
            conversion_table['T'] = conversion_table['t'] = 3;
            conversion_table['U'] = conversion_table['T']; // set U equal to T
            conversion_table['u'] = conversion_table['t']; // set u equal to t

            return conversion_table;
        }()
    };

    // The forward declaration of the complement table. With this, my_dna4 models seqan3::nucleotide_alphabet.
    static const std::array<my_dna4, alphabet_size> complement_table;

    friend nucleotide_base<my_dna4, 4>; // Grant seqan3::nucleotide_base access to private/protected members.
    friend nucleotide_base<my_dna4, 4>::base_t; // Grant seqan3::alphabet_base access to private/protected members.
};

// Defines the `_my_dna4` *char literal* so you can write `'C'_my_dna4` instead of `my_dna4{}.assign_char('C')`.
constexpr my_dna4 operator""_my_dna4(char const c) noexcept
{
    return my_dna4{}.assign_char(c);
}

// The definition of the complement table. With this, my_dna4 models seqan3::nucleotide_alphabet.
constexpr std::array<my_dna4, my_dna4::alphabet_size> my_dna4::complement_table
{
    'T'_my_dna4,    // the complement of 'A'_my_dna4
    'G'_my_dna4,    // the complement of 'C'_my_dna4
    'C'_my_dna4,    // the complement of 'G'_my_dna4
    'A'_my_dna4     // the complement of 'T'_my_dna4
};

int main()
{
    my_dna4 my_letter{'C'_my_dna4};

    my_letter.assign_char('S'); // Characters other than A,C,G,T are implicitly converted to `A`.

    seqan3::debug_stream << my_letter; // "A";
}
