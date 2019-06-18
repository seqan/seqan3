#include <seqan3/core/debug_stream.hpp>
#include <seqan3/std/span>

using namespace seqan3;

int main()
{
    std::string text{"Garfield the fat cat without a hat."};
    size_t start{2};
    size_t span{3};

    std::span text_view{std::data(text) + start, span}; // represent interval [2, 4]

    debug_stream << text_view << '\n'; // Prints "rfi"
}
