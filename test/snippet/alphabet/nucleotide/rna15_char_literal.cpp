#include <seqan3/alphabet/nucleotide/rna15.hpp>

int main()
{
    using namespace seqan3::literals;

    seqan3::rna15 letter1{'A'_rna15};
    auto letter2 = 'A'_rna15;
}
