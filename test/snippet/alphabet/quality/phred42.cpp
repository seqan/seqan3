#include <seqan3/alphabet/nucleotide/dna4.hpp>
#include <seqan3/alphabet/quality/phred42.hpp>
#include <seqan3/alphabet/quality/qualified.hpp>
#include <seqan3/core/debug_stream.hpp>

using seqan3::operator""_dna4;
using seqan3::operator""_phred42;

int main()
{
    seqan3::phred42 phred;
    phred.assign_phred(2);
    seqan3::debug_stream << phred.to_phred() << "\n"; // 2
    seqan3::debug_stream << phred.to_char() << "\n";  // '#'
    seqan3::debug_stream << phred.to_rank() << "\n";  // 2

    seqan3::phred42 another_phred = '('_phred42;
    seqan3::debug_stream << another_phred.to_phred() << "\n"; // 52

    another_phred.assign_phred(49); // converted down to 41
    seqan3::debug_stream << another_phred.to_phred() << "\n"; // 41

    std::vector<seqan3::qualified<seqan3::dna4, seqan3::phred42>> query{{'A'_dna4, '!'_phred42},
                                                                        {'C'_dna4, 'A'_phred42},
                                                                        {'G'_dna4, '6'_phred42},
                                                                        {'T'_dna4, '&'_phred42}};
}
