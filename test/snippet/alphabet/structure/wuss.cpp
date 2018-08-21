#include <seqan3/alphabet/structure/wuss.hpp>

int main()
{
using namespace seqan3;
{
//! [general]
// create vector
std::vector<wuss51> vec{wuss51::UNPAIRED, wuss51::PAIR_CLOSE, wuss51::PAIR_CLOSE};
// modify and print
vec[1] = wuss51::PAIR_OPEN;
for (wuss51 chr : vec)
    std::cout << chr;  // .<>
std::cout << "\n";
//! [general]
}
}
