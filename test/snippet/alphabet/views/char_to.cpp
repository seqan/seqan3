#include <seqan3/alphabet/nucleotide/dna4.hpp>
#include <seqan3/alphabet/nucleotide/dna5.hpp>
#include <seqan3/alphabet/views/char_to.hpp>
#include <seqan3/core/debug_stream.hpp>

int main()
{
    std::string str{"ACTTTGATAN"};
    seqan3::debug_stream << (str | seqan3::views::char_to<seqan3::dna4>) << '\n'; // ACTTTGATAA
    seqan3::debug_stream << (str | seqan3::views::char_to<seqan3::dna5>) << '\n'; // ACTTTGATAN
}
