#include <seqan3/test/snippet/create_temporary_snippet_file.hpp>
create_temporary_snippet_file my_fasta
{
    "my.fasta",
R"//![fasta_file](
>seq1
AVAV
>seq2
AVAVA
)//![fasta_file]"
}; // std::filesystem::current_path() / "my.fasta" will be deleted after the execution

//![main]
#include <seqan3/io/sequence_file/all.hpp>

int main()
{
    using traits_type = seqan3::sequence_file_input_default_traits_aa;
    seqan3::sequence_file_input<traits_type> fin{std::filesystem::current_path() / "my.fasta"};
}
//![main]
