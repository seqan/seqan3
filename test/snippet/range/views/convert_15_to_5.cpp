#include <seqan3/alphabet/nucleotide/dna15.hpp>
#include <seqan3/alphabet/nucleotide/dna5.hpp>
#include <seqan3/range/views/convert.hpp>

int main()
{
    using seqan3::operator""_dna15;

    seqan3::dna15_vector vec2{"ACYGTN"_dna15};
    auto v4 = vec2 | seqan3::views::convert<seqan3::dna5>; // == "ACNGTN"_dna5
}
