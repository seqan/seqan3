#include <iostream>                     // for std::cout
#include <seqan3/alphabet/all.hpp>      // include all alphabet headers

template <seqan3::alphabet t>
void print(t const v)
{
    std::cout << "I am an alphabet and my value as char is: " << seqan3::to_char(v) << '\n';
}

using seqan3::operator""_dna5;
using seqan3::operator""_aa27;

int main()
{
    auto d = 'A'_dna5;
    auto a = 'L'_aa27;
    auto g = seqan3::gap{};

    print(d);
    print(a);
    print(g);
}
