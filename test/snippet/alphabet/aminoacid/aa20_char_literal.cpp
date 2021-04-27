#include <seqan3/alphabet/aminoacid/aa20.hpp>

int main()
{
    using namespace seqan3::literals;

    seqan3::aa20 letter1{'A'_aa20};
    auto letter2 = 'A'_aa20;
}
