#include <seqan3/test/snippet/create_temporary_snippet_file.hpp>
create_temporary_snippet_file my_fastq
{
    "my.fastq",
R"//![fastq_file](
@seq1
CGATCGATC
+
IIIIIIIII
@seq2
AGCG
+
IIII
@seq3
AGCTAGCAGCGATCG
+
IIIIIHIIJJIIIII
@seq4
AGC
+
III
@seq5
AGCTAGCAGCGATCG
+
IIIIIHIIJJIIIII
)//![fastq_file]"
}; // std::filesystem::current_path() / "my.fastq" will be deleted after the execution

// std::filesystem::current_path() / "output.fastq" will be deleted after the execution
create_temporary_snippet_file output{"output.fastq", ""};

//![solution]
#include <seqan3/io/sequence_file/all.hpp>
#include <seqan3/std/filesystem>
#include <seqan3/std/ranges>

int main()
{
    std::filesystem::path current_path = std::filesystem::current_path();

    seqan3::sequence_file_input fin{current_path / "my.fastq"};
    seqan3::sequence_file_output fout{current_path / "output.fastq"};

    auto length_filter = std::views::filter([] (auto const & rec)
    {
        return std::ranges::size(rec.sequence()) >= 5;
    });

    for (auto & record : fin | length_filter)
    {
        fout.push_back(record);
    }
}
//![solution]
