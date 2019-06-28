#include <seqan3/alphabet/nucleotide/dna4.hpp>
#include <seqan3/core/debug_stream.hpp>

int main()
{
    seqan3::dna4 my_letter;
    seqan3::assign_rank_to(0, my_letter);               // assign an A via rank interface
    seqan3::assign_char_to('A', my_letter);             // assign an A via char interface

    std::cout << seqan3::to_char(my_letter);            // prints 'A'
    std::cout << (unsigned)seqan3::to_rank(my_letter);  // prints 0
    // we have to add the cast here, because uint8_t is also treated as a char type by default :(

    // Using SeqAn's debug_stream:
    seqan3::debug_stream << seqan3::to_char(my_letter); // prints 'A'
    seqan3::debug_stream << my_letter;                  // prints 'A' (calls to_char() automatically!)
    seqan3::debug_stream << seqan3::to_rank(my_letter); // prints 0   (casts uint8_t to unsigned automatically!)
}
