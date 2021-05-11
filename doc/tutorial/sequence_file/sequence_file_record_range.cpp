#include <seqan3/test/snippet/create_temporary_snippet_file.hpp>
// std::filesystem::current_path() / "my.fastq" will be deleted after the execution
create_temporary_snippet_file my_fastq{"my.fasta", "\n"};

//![include]
#include <seqan3/io/sequence_file/all.hpp>
//![include]

int main()
{
    seqan3::sequence_file_input file{std::filesystem::current_path() / "my.fasta"};

//![record_range]
for (auto & record : file)
{
    // do something with my record
}
//![record_range]
}
