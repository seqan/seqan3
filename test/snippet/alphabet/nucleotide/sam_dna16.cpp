#include <seqan3/alphabet/nucleotide/sam_dna16.hpp>
#include <seqan3/core/detail/reflection.hpp>
#include <seqan3/core/debug_stream.hpp>

using namespace seqan3;

int main()
{
    sam_dna16 my_letter{'A'_sam_dna16};
    // doesn't work:
    // sam_dna16 my_letter{'A'};

    my_letter.assign_char('='); // <- this does!

    my_letter.assign_char('F'); // unknown characters are implicitly converted to N.
    debug_stream << my_letter << std::endl; // "N";

    //! [operator""_sam_dna16]
    // these don't work:
    // sam_dna16_vector foo{"ACGTTA"};
    // sam_dna16_vector bar = "ACGTTA";

    // but these do:
    sam_dna16_vector foo{"ACgtTA"_sam_dna16};
    sam_dna16_vector bar = "ACG==A"_sam_dna16;
    auto bax = "A=GTT!"_sam_dna16;

    debug_stream << foo << "\n" << bar << "\n" << bax << "\n";
    //! [operator""_sam_dna16]
}
