#include <seqan3/alphabet/aminoacid/aa10li.hpp>

int main()
{
    using namespace seqan3::literals;

    seqan3::aa10li letter1{'A'_aa10li};
    auto letter2 = 'A'_aa10li;
}
