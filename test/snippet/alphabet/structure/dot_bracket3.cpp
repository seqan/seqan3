#include <seqan3/alphabet/structure/dot_bracket3.hpp>

using namespace seqan3;

int main()
{
//! [general]
// create vector
std::vector<dot_bracket3> vec{dot_bracket3::UNPAIRED, dot_bracket3::PAIR_CLOSE, dot_bracket3::PAIR_CLOSE};
// modify and print
vec[1] = dot_bracket3::PAIR_OPEN;
for (dot_bracket3 chr : vec)
    std::cout << chr;  // .()
std::cout << "\n";
//! [general]
}
