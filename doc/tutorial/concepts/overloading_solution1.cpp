#include <iostream>                     // for std::cout
#include <seqan3/alphabet/all.hpp>      // include all alphabet headers

using namespace seqan3;

template <Alphabet t>
void print(t const v)
{
    std::cout << "I am an alphabet and my value as char is: " << to_char(v) << '\n';
}

int main()
{
    auto d = 'A'_dna5;
    auto a = 'L'_aa27;
    auto g = gap{};

    print(d);
    print(a);
    print(g);
}
