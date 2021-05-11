#include <seqan3/test/snippet/create_temporary_snippet_file.hpp>
create_temporary_snippet_file my_fastq
{
    "my.fastq",
R"//![fastq_file](
@seq1
AGCTAGCAGCGATCG
+
IIIIIHIIIIIIIII
@seq2
CGATCGATC
+
IIIIIIIII
@seq3
AGCGATCGAGGAATATAT
+
IIIIHHGIIIIHHGIIIH
)//![fastq_file]"
}; // std::filesystem::current_path() / "my.fastq" will be deleted after the execution

// std::filesystem::current_path() / "output.fasta" will be deleted after the execution
create_temporary_snippet_file output_fasta{"output.fasta", ""};

//![main]
#include <seqan3/io/sequence_file/all.hpp>

int main()
{
    auto current_path = std::filesystem::current_path();

    seqan3::sequence_file_output{current_path / "output.fasta"}
        = seqan3::sequence_file_input{current_path / "my.fastq"};
}
//![main]
