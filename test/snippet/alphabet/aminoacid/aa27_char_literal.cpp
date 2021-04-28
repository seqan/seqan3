#include <seqan3/alphabet/aminoacid/aa27.hpp>

int main()
{
    using namespace seqan3::literals;

    seqan3::aa27 letter1{'A'_aa27};
    auto letter2 = 'A'_aa27;
}
