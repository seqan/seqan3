#include <iostream>

#include <seqan3/range/view/single_pass_input.hpp>

using namespace seqan3;

int main()
{
//! [usage]
std::string str{"hello"};
auto v = str | view::single_pass_input;
std::cout << *++v.begin() << std::endl;  // prints 'e'
std::cout << *++v.begin() << std::endl;  // prints 'l'
std::cout << *++v.begin() << std::endl;  // prints 'l'
std::cout << *++v.begin() << std::endl;  // prints 'o'
//! [usage]
}
