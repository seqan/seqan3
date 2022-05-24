#include <seqan3/alphabet/nucleotide/dna4.hpp>
#include <seqan3/alphabet/views/char_strictly_to.hpp>
#include <seqan3/core/debug_stream.hpp>

int main()
{
    std::string_view str{"ACTTTGATAN"};
    try
    {
        seqan3::debug_stream << (str | seqan3::views::char_strictly_to<seqan3::dna4>); // ACTTTGATA
    }
    catch (seqan3::invalid_char_assignment const &)
    {
        seqan3::debug_stream << "\n[ERROR] Invalid char!\n"; // Will throw on parsing 'N'
    }
}
