#include <seqan3/alphabet/structure/dssp9.hpp>

int main()
{
    using namespace seqan3::literals;

    seqan3::dssp9 letter1{'('_dssp9};
    auto letter2 = '('_dssp9;
}
