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

//![main]
#include <seqan3/io/sequence_file/all.hpp>
#include <seqan3/utility/views/chunk.hpp>

int main()
{
    seqan3::sequence_file_input fin{std::filesystem::current_path() / "my.fastq"};

    // `&&` is important because seqan3::views::chunk returns temporaries!
    for (auto && records : fin | seqan3::views::chunk(10))
    {
        // `records` contains 10 elements (or less at the end)
        seqan3::debug_stream << "Taking the next 10 sequences:\n";
        seqan3::debug_stream << "ID:  " << (*records.begin()).id() << '\n'; // prints first ID in batch
    }
//![main]
}
