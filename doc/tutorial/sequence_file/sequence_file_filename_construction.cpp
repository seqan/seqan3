// this macro will rename the `main` function below as `sam_file_filename_construction`.
// we just want to show that the syntax does compile.
#define main sequence_file_filename_construction

//![main]
#include <seqan3/std/filesystem>

#include <seqan3/io/sequence_file/all.hpp>

int main()
{
    auto filename = std::filesystem::current_path() / "my.fasta";

    seqan3::sequence_file_input fin_from_filename{filename};

    seqan3::sequence_file_input fin_from_stream{std::cin, seqan3::format_fasta{}};

    return 0;
}
//![main]

// this makes the snippet an executable
#undef main
int main() {} // do nothing
