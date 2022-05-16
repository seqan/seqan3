#include <iostream>                     // for std::cout
#include <seqan3/alphabet/all.hpp>      // include all alphabet headers

template <seqan3::alphabet t>
void print(t const v)
{
    std::cout << "I am an alphabet and my value as char is: " << seqan3::to_char(v) << '\n';
}


int main()
{
    using namespace seqan3::literals;

    auto d = 'A'_dna5;
    auto a = 'L'_aa27;
    auto g = seqan3::gap{};

    print(d);
    print(a);
    print(g);
}
