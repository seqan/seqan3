#include <seqan3/core/debug_stream.hpp>       // pretty printing
#include <seqan3/search/algorithm/search.hpp> // this already includes (bi)_fm_index

using namespace seqan3;
using namespace std::string_literals; // for using the ""s string literal

int main()
{
    std::string text{"Garfield the fat cat without a hat."};
    fm_index index{text};
    debug_stream << search("cat"s, index) << '\n'; // [17]
}
