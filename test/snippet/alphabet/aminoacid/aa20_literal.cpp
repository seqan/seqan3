#include <seqan3/alphabet/aminoacid/aa20.hpp>

int main()
{
    using seqan3::operator""_aa20;
    
    seqan3::aa20_vector foo{"ABFUYR"_aa20};
    seqan3::aa20_vector bar = "ABFUYR"_aa20;
    auto bax = "ABFUYR"_aa20;
}
