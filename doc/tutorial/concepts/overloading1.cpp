#include <iostream>                 // for std::cout

#include <seqan3/std/concepts>      // GCC7 - GCC9 or
//#include <concepts>               // compilers with full C++20 support

template <std::Integral t>
void print(t const v)
{
    std::cout << "Integral value: " << v << '\n';
}

int main()
{
    int i{4};
    unsigned u{3};

    print(i);                       // prints "Integral value: 4"
    print(u);                       // prints "Integral value: 3"
}
