#include <vector>

#include <seqan3/alphabet/structure/dssp9.hpp>
#include <seqan3/core/debug_stream.hpp>

int main()
{
    using seqan3::operator""_dssp9;

    // create vector
    std::vector<seqan3::dssp9> vec{'E'_dssp9, 'H'_dssp9, 'H'_dssp9, 'H'_dssp9, 'T'_dssp9, 'G'_dssp9};

    // modify and print
    vec[1] = 'C'_dssp9;

    for (seqan3::dssp9 chr : vec)
        seqan3::debug_stream << seqan3::to_char(chr);  // ECHHTG
        
    seqan3::debug_stream << "\n";
}
