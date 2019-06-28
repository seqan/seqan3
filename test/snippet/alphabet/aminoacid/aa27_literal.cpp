#include <seqan3/alphabet/aminoacid/aa27.hpp>

int main()
{
    using seqan3::operator""_aa27;
    
    seqan3::aa27_vector foo{"ABFUYR"_aa27};
    seqan3::aa27_vector bar = "ABFUYR"_aa27;
    auto bax = "ABFUYR"_aa27;
}
