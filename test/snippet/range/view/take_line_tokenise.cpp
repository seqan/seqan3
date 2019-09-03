#include <seqan3/core/debug_stream.hpp>
#include <seqan3/range/view/take_line.hpp>
#include <seqan3/range/view/single_pass_input.hpp>

int main()
{
    std::string vec{"foo\nbar"};
    auto v = vec | seqan3::view::single_pass_input | seqan3::view::take_line;
    seqan3::debug_stream << v << '\n'; // [f,o,o]
    seqan3::debug_stream << v << '\n'; // [b,a,r] (parsing it again gives us the next line)
}
