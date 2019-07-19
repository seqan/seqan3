#include <seqan3/alphabet/nucleotide/dna3bs.hpp>
#include <seqan3/core/debug_stream.hpp>

int main()
{
    using seqan3::operator""_dna3bs;

    seqan3::dna3bs_vector foo{"ACGTTA"_dna3bs}; // ATGTTA
    seqan3::dna3bs_vector bar = "ATGTTA"_dna3bs;

    if (foo == bar)
        seqan3::debug_stream << "yeah\n"; // "yeah";

    auto bax = "ACGTTA"_dna3bs;

    return 0;
}
