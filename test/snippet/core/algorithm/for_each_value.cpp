#include <iostream>
#include <string>

#include <seqan3/core/algorithm/parameter_pack.hpp>
#include <seqan3/core/debug_stream.hpp>

using namespace seqan3;

int main()
{
    auto fn = [](auto && a)
    {
        debug_stream << a;
    };

    // prints each argument, i.e. "0, 1, 2, 3\n"
    detail::for_each_value(fn, 0, ", ", 1.0, ", ", std::string{"2, 3"}, '\n');

    // is the same as explicitly writing
    fn(0);
    fn(", ");
    fn(1.0);
    fn(", ");
    fn(std::string{"2, 3"});
    fn('\n');
    return 0;
}
