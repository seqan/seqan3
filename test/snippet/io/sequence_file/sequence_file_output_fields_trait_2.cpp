#include <sstream>

#include <seqan3/io/sequence_file/input.hpp>
#include <seqan3/io/sequence_file/output.hpp>

auto input = R"(@TEST1
ACGT
+
##!#
@Test2
AGGCTGA
+
##!#!!!
@Test3
GGAGTATAATATATATATATATAT
+
##!###!###!###!###!###!#)";

int main()
{
    seqan3::sequence_file_input  fin{std::istringstream{input},
                                     seqan3::format_fastq{},
                                     seqan3::fields<seqan3::field::SEQ, seqan3::field::ID, seqan3::field::QUAL>{}};
    seqan3::sequence_file_output fout{std::ostringstream{},
                                      seqan3::format_fastq{}}; // doesn't have to match the configuration

    for (auto & r : fin)
    {
        if (true) // r fulfills some criterium
            fout.push_back(r);
    }
}
