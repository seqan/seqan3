#include <seqan3/alphabet/nucleotide/sam_dna16.hpp>
#include <seqan3/core/debug_stream.hpp>

int main()
{
    using seqan3::operator""_sam_dna16;
    
    seqan3::sam_dna16_vector foo{"ACgtTA"_sam_dna16};
    seqan3::sam_dna16_vector bar = "ACG==A"_sam_dna16;
    auto bax = "A=GTT!"_sam_dna16;

    seqan3::debug_stream << foo << "\n" << bar << "\n" << bax << "\n";
}
