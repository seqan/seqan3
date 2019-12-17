#include <sstream>
#include <string>
#include <tuple>
#include <vector>

#include <seqan3/alphabet/nucleotide/dna5.hpp>
#include <seqan3/alphabet/quality/phred42.hpp>
#include <seqan3/alphabet/quality/qualified.hpp>
#include <seqan3/io/sequence_file/output.hpp>

int main()
{
    using seqan3::operator""_dna5;
    using seqan3::operator""_phred42;

    seqan3::sequence_file_output fout{std::ostringstream{},
                                      seqan3::format_fasta{},
                                      seqan3::fields<seqan3::field::id, seqan3::field::seq_qual>{}};

    for (int i = 0; i < 5; i++)
    {
        std::string id{"test_id"};
        // vector of combined data structure:
        std::vector<seqan3::qualified<seqan3::dna5, seqan3::phred42>> seq_qual{{'N'_dna5, '7'_phred42}};

        // ...

        fout.emplace_back(id, seq_qual);       // note also that the order the argumets is now different, because
        // or:                                    you specified that ID should be first in the fields template argument
        fout.push_back(std::tie(id, seq_qual));
    }
}
