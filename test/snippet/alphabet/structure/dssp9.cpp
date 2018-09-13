#include <seqan3/alphabet/structure/dssp9.hpp>

using namespace seqan3;

int main()
{
//! [general]
// create vector
std::vector<dssp9> vec{dssp9::E, dssp9::H, dssp9::H, dssp9::H, dssp9::T, dssp9::G};
// modify and print
vec[1] = dssp9::C;
for (dssp9 chr : vec)
    std::cout << to_char(chr);  // ECHHTG
std::cout << "\n";
//! [general]
}
