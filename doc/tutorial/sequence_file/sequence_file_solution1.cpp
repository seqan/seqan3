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

//![solution]
#include <seqan3/std/filesystem>

//![include_debug_stream]
#include <seqan3/core/debug_stream.hpp>
//![include_debug_stream]
#include <seqan3/io/sequence_file/all.hpp>

int main()
{
    std::filesystem::path current_path = std::filesystem::current_path();

    seqan3::sequence_file_input fin{current_path / "my.fastq"};

    for (auto & rec : fin)
    {
        seqan3::debug_stream << "ID:  "  << rec.id() << '\n';
        seqan3::debug_stream << "SEQ: "  << rec.sequence() << '\n';
        seqan3::debug_stream << "QUAL: " << rec.base_qualities() << '\n';
    }
}
//![solution]
