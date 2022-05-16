#include <concepts>
#include <iostream> // for std::cout

template <std::integral t>
void print(t const v)
{
    std::cout << "integral value: " << v << '\n';
}

int main()
{
    int i{4};
    unsigned u{3};

    print(i); // prints "integral value: 4"
    print(u); // prints "integral value: 3"
}
