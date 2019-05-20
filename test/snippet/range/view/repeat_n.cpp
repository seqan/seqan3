#include <seqan3/core/debug_stream.hpp>
#include <seqan3/range/view/repeat_n.hpp>
#include <seqan3/std/ranges>

using namespace seqan3;

int main()
{
    auto v = view::repeat_n(std::string{"foo"}, 5);

    debug_stream << v.size() << std::endl; // prints 5
    debug_stream << v << std::endl;        // prints ["foo", "foo", "foo", "foo", "foo"]

    v[0] = std::string{"foobar"};

    // Now it always prints "foobar"
    debug_stream << *std::ranges::begin(v) << std::endl; // prints "foobar"
    debug_stream << v.size() << std::endl;               // prints 5; Note: the size cannot be changed anymore
}
