#include <seqan3/alphabet/structure/wuss.hpp>
#include <seqan3/io/stream/debug_stream.hpp>

using namespace seqan3;

int main()
{
//! [general]
// create vector
std::vector<wuss51> vec{wuss51::UNPAIRED, wuss51::PAIR_CLOSE, wuss51::PAIR_CLOSE};
// modify and print
vec[1] = wuss51::PAIR_OPEN;
for (wuss51 chr : vec)
    debug_stream << to_char(chr);  // .<>
debug_stream << "\n";
//! [general]
}
