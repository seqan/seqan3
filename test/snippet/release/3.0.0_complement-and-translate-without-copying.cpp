#include <vector>

#include <seqan3/alphabet/nucleotide/dna5.hpp>
#include <seqan3/alphabet/views/complement.hpp>
#include <seqan3/alphabet/views/translate.hpp>
#include <seqan3/core/debug_stream.hpp>

int main()
{
    using namespace seqan3::literals;

    seqan3::dna5_vector sequence{"ACGTACGTACGTA"_dna5};
    // pipe the vector through two view adaptors:
    auto aminoacid_sequence = sequence | seqan3::views::complement | seqan3::views::translate_single;
    // aminoacid_sequence is a view of length 4, accessing the elements on-demand will return

    seqan3::debug_stream << aminoacid_sequence << '\n'; // "CMHA"
}
