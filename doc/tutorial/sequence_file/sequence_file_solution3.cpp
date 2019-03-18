#include <fstream>

#include <seqan3/std/filesystem>

struct write_file_dummy_struct
{
    write_file_dummy_struct()
    {

auto file_raw = R"//![fastq_file](
@seq1
CGATCGATC
+
IIIIIIIII
@seq2
AGCGATCGAGGAATATAT
+
IIIIHHGIIIIHHGIIIH
@seq3
AGCTAGCAGCGATCG
+
IIIIIHIIJJIIIII
@seq4
AGCGATCGAGGAATATAT
+
IIIIHHGIIIIHHGIIIH
@seq5
AGCTAGCAGCGATCG
+
IIIIIHIIJJIIIII
)//![fastq_file]";

        std::ofstream file{std::filesystem::temp_directory_path()/"my.fastq"};
        std::string str{file_raw};
        file << str.substr(1); // skip first newline
    }
};

write_file_dummy_struct go{};

//![solution]
#include <range/v3/numeric/accumulate.hpp>   // ranges::accumulate
#include <range/v3/view/take.hpp>            // ranges::take

#include <seqan3/io/sequence_file/all.hpp>
#include <seqan3/io/stream/debug_stream.hpp>
#include <seqan3/std/ranges>
#include <seqan3/std/filesystem>

using namespace seqan3;

int main()
{
    std::filesystem::path tmp_dir = std::filesystem::temp_directory_path(); // get the temp directory

    sequence_file_input fin{tmp_dir/"my.fastq"};

    auto minimum_quality_filter = std::view::filter([] (auto const & rec)
    {
        auto qual = get<field::QUAL>(rec) | std::view::transform([] (auto q) { return q.to_phred(); });
        double sum = ranges::accumulate(qual.begin(), qual.end(), 0);
        return sum / std::ranges::size(qual) >= 40;
    });

    for (auto & rec : fin | minimum_quality_filter | std::view::take(2))
    {
        debug_stream << "ID: " << get<field::ID>(rec) << '\n';
    }
}
//![solution]
