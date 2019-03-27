#include <vector>

#include <seqan3/io/stream/debug_stream.hpp>
#include <seqan3/range/view/to_lower.hpp>
#include <seqan3/range/view/to_upper.hpp>

using namespace seqan3;

int main()
{

{
//! [to_lower]
std::string s{"CHanGED!"};
std::vector<std::string> sv{"CHANGED", "unchanged!"};
auto v1 = s | view::to_lower; 
auto v2 = sv | view::to_lower; 

debug_stream << v1 << '\n'; // => "changed!"
debug_stream << v2 << '\n'; // => ["changed", "unchanged!"]
//! [to_lower]

}

{
//! [to_upper]
std::string s{"CHanGED!"};
std::vector<std::string> sv{"changed", "UNCHANGED!"};
auto v1 = s | view::to_upper; 
auto v2 = sv | view::to_upper; 

debug_stream << v1 << '\n'; // => "CHANGED!"
debug_stream << v2 << '\n'; // => ["CHANGED", "UNCHANGED!"]
//! [to_upper]
}

}
