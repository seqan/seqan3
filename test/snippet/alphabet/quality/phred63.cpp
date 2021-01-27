#include <seqan3/alphabet/nucleotide/dna4.hpp>
#include <seqan3/alphabet/quality/phred63.hpp>
#include <seqan3/alphabet/quality/qualified.hpp>
#include <seqan3/core/debug_stream.hpp>

using seqan3::operator""_dna4;
using seqan3::operator""_phred63;

int main()
{
    seqan3::phred63 phred;
    phred.assign_rank(2); // wrapper for assign_phred(2)
    seqan3::debug_stream << phred.to_phred() << "\n"; // 2
    seqan3::debug_stream << phred.to_char() << "\n";        // '#'
    seqan3::debug_stream << phred.to_rank() << "\n";  // 2

    seqan3::phred63 another_phred{49};
    seqan3::debug_stream << another_phred.to_phred() << "\n"; // 49

    seqan3::phred63 a_third_phred;
    another_phred.assign_phred(75); // converted down to 62
    seqan3::debug_stream << another_phred.to_phred() << "\n"; // 62

    std::vector<seqan3::qualified<seqan3::dna4, seqan3::phred63>> query{{'A'_dna4, '!'_phred63},
                                                                        {'C'_dna4, 'A'_phred63},
                                                                        {'G'_dna4, '6'_phred63},
                                                                        {'T'_dna4, '&'_phred63}};
}
