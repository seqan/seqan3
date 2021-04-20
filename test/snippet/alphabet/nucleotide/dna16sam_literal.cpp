#include <seqan3/alphabet/nucleotide/dna16sam.hpp>
#include <seqan3/core/debug_stream.hpp>

int main()
{
    using seqan3::operator""_dna16sam;

    seqan3::dna16sam_vector foo{"ACgtTA"_dna16sam};
    seqan3::dna16sam_vector bar = "ACG==A"_dna16sam;
    auto bax = "A=GTT!"_dna16sam;

    seqan3::debug_stream << foo << "\n" << bar << "\n" << bax << "\n";
}
