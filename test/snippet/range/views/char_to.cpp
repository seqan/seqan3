#include <string>

#include <seqan3/alphabet/nucleotide/dna4.hpp>
#include <seqan3/alphabet/nucleotide/dna5.hpp>
#include <seqan3/alphabet/views/char_to.hpp>

int main()
{
    std::string s{"ACTTTGATAN"};
    auto v1 = s | seqan3::views::char_to<seqan3::dna4>; // == "ACTTTGATAA"_dna4
    auto v2 = s | seqan3::views::char_to<seqan3::dna5>; // == "ACTTTGATAN"_dna5
}
