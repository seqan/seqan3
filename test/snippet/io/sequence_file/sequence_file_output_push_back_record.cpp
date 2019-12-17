#include <sstream>
#include <string>

#include <seqan3/alphabet/nucleotide/dna5.hpp>
#include <seqan3/core/type_list/type_list.hpp>
#include <seqan3/io/sequence_file/output.hpp>

int main()
{
    using seqan3::operator""_dna5;

    seqan3::sequence_file_output fout{std::ostringstream{}, seqan3::format_fasta{}};
    for(int i = 0; i < 5; ++i) // some criteria
    {
        seqan3::record<seqan3::type_list<seqan3::dna5_vector, std::string>,
                       seqan3::fields<seqan3::field::seq, seqan3::field::id>> r{"ACGT"_dna5, "ID1"};

        // ...

        fout.push_back(r);
    }
}
