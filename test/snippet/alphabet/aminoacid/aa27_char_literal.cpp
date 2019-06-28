#include <seqan3/alphabet/aminoacid/aa27.hpp>

int main()
{
    using seqan3::operator""_aa27;

    seqan3::aa27 acid1{'A'_aa27};
    
    auto acid2 = 'Y'_aa27; // type = aa27
}
