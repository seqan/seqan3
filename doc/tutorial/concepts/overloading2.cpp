#include <iostream>                 // for std::cout

#include <seqan3/std/concepts>      // GCC7 - GCC9 or
//#include <concepts>               // compilers with full C++20 support

template <std::integral t>
void print(t const v)
{
    std::cout << "integral value: " << v << '\n';
}

template <std::unsigned_integral t>
void print(t const v)
{
    std::cout << "Unsigned value: " << v << '\n';
}

int main()
{
    int i{4};
    unsigned u{3};

    print(i);                       // prints "integral value: 4"
    print(u);                       // prints "Unsigned value: 3"
}
