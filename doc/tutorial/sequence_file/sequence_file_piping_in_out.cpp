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

// std::filesystem::current_path() / "output.fastq" will be deleted after the execution
create_temporary_snippet_file output_fastq{"output.fastq", ""};

//![main]
#include <seqan3/std/algorithm>
#include <seqan3/std/ranges>

#include <seqan3/io/sequence_file/all.hpp>

int main()
{
    auto current_path = std::filesystem::current_path();

    seqan3::sequence_file_input fin{current_path / "my.fastq"};
    seqan3::sequence_file_output fout{current_path / "output.fastq"};

    // the following are equivalent:
    // 1. copy records of input file into output file
    std::ranges::move(fin, fout.begin());

    // 2. assign all records of input file to output file
    fout = fin;

    // 3. same as 2. but as one liner
    seqan3::sequence_file_output{current_path / "output.fastq"}
        = seqan3::sequence_file_input{current_path / "my.fastq"};
}
//![main]
