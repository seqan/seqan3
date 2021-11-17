#include <seqan3/test/snippet/create_temporary_snippet_file.hpp>
// std::filesystem::current_path() / "my.fasta" will be deleted after the execution
seqan3::test::create_temporary_snippet_file my_fasta{"my.fasta", ""};

//![main]
#include <filesystem>

#include <seqan3/io/sequence_file/output.hpp>

int main()
{
    auto fasta_file = std::filesystem::current_path() / "my.fasta";

    // FASTA format detected, std::ofstream opened for file
    seqan3::sequence_file_output fin{fasta_file};
}
//![main]
