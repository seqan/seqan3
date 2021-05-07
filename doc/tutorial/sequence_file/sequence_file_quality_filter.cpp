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
#include <numeric> // std::accumulate
#include <seqan3/std/ranges>

#include <seqan3/core/debug_stream.hpp>
#include <seqan3/io/sequence_file/all.hpp>

int main()
{
    seqan3::sequence_file_input fin{std::filesystem::current_path() / "my.fastq"};

    // std::views::filter takes a function object (a lambda in this case) as input that returns a boolean
    auto minimum_quality_filter = std::views::filter([] (auto const & rec)
    {
        auto qualities = rec.base_qualities()
                       | std::views::transform([] (auto quality) { return seqan3::to_phred(quality); });

        auto sum = std::accumulate(qualities.begin(), qualities.end(), 0);
        return sum / std::ranges::size(qualities) >= 40; // minimum average quality >= 40
    });

    for (auto & rec : fin | minimum_quality_filter)
    {
        seqan3::debug_stream << "ID: " << rec.id() << '\n';
    }
}
//![main]
