#include <iostream>

#include <seqan3/core/debug_stream.hpp>
#include <seqan3/range/view/single_pass_input.hpp>

using namespace seqan3;

int main()
{
//! [usage]
std::string str{"hello"};
auto v = str | view::single_pass_input;
debug_stream << *++v.begin() << std::endl;  // prints 'e'
debug_stream << *++v.begin() << std::endl;  // prints 'l'
debug_stream << *++v.begin() << std::endl;  // prints 'l'
debug_stream << *++v.begin() << std::endl;  // prints 'o'
//! [usage]
}
