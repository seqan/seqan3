#include <string>

#include <seqan3/core/debug_stream.hpp>
#include <seqan3/std/ranges>

int main()
{
//![usage]
std::views::take_while([] (auto const & l) { return (l != '\r') && (l != '\n'); });
//![usage]
}
