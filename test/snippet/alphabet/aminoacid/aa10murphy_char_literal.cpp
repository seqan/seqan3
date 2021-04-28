#include <seqan3/alphabet/aminoacid/aa10murphy.hpp>

int main()
{
    using namespace seqan3::literals;

    seqan3::aa10murphy letter1{'A'_aa10murphy};
    auto letter2 = 'A'_aa10murphy;
}
