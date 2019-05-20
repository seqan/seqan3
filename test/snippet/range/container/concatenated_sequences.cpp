#include <seqan3/alphabet/nucleotide/dna4.hpp>
#include <seqan3/core/debug_stream.hpp>
#include <seqan3/range/container/concatenated_sequences.hpp>

using namespace seqan3;

int main()
{

{
std::vector<dna4> vector_of_length1000;
bool not_full = false;
//! [usage]
concatenated_sequences<dna4_vector> concat1{"ACGT"_dna4, "GAGGA"_dna4};
debug_stream << concat1[0] << '\n'; // "ACGT"

std::vector<dna4_vector> concat2{"ACTA"_dna4, "AGGA"_dna4};

concat1 = concat2;               // you can assign from other ranges

concat2[0] = "ATTA"_dna4;        // this works for vector of vector
//concat1[0] = "ATTA"_dna4;      // but not on concatenated_sequences

concat1[0][1] = 'T'_dna4;         // this, however, does
debug_stream << concat1[0] << '\n'; // "ATTA"


// if you know that you will be adding a thousand vectors of length thousand:
concat1.reserve(1'000);
concat1.concat_reserve(1'000 * 1'000);
while (not_full)
{
    // ...
    concat1.push_back(vector_of_length1000);
}
//! [usage]
}

{
//! [insert]
concatenated_sequences<dna4_vector> foobar;
foobar.insert(foobar.end(), "ACGT"_dna4);
debug_stream << foobar[0] << '\n'; // "ACGT"
//! [insert]
}

{
//! [insert2]
concatenated_sequences<dna4_vector> foobar;
foobar.insert(foobar.end(), 2, "ACGT"_dna4);
debug_stream << foobar[0] << '\n'; // "ACGT"
debug_stream << foobar[1] << '\n'; // "ACGT"
//! [insert2]
}
}
