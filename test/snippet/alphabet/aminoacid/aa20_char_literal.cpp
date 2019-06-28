#include <seqan3/alphabet/aminoacid/aa20.hpp>

int main()
{
    using seqan3::operator""_aa20;

    seqan3::aa20 acid1{'A'_aa20};
    
    auto acid2 = 'Y'_aa20; // type = aa20
}
