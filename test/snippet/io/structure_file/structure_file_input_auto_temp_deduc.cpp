#include <sstream>

#include <seqan3/alphabet/nucleotide/rna4.hpp>
#include <seqan3/alphabet/structure/wuss.hpp>
#include <seqan3/io/structure_file/input.hpp>
#include <seqan3/io/structure_file/output.hpp>
#include <seqan3/std/filesystem>

auto input = R"(> S.cerevisiae_tRNA-PHE M10740/1-73
GCGGAUUUAGCUCAGUUGGGAGAGCGCCAGACUGAAGAUUUGGAGGUCCUGUGUUCGAUCCACAGAAUUCGCA
(((((((..((((........)))).((((.........)))).....(((((.......)))))))))))). (-17.50)
> example
UUGGAGUACACAACCUGUACACUCUUUC
..(((((..(((...)))..)))))... (-3.71))";

int main()
{
    auto tmp_file = std::filesystem::temp_directory_path() / "my.dbn";

    using namespace seqan3::literals;

    // First, make /tmp/input.dbn
    {
        seqan3::structure_file_output fout{tmp_file};
        fout.emplace_back("GCGGAUUUAGCUCAGUUGGGAGAGCGCCAGACUGAAGAUUUGGAGGUCCUGUGUUCGAUCCACAGAAUUCGCA"_rna4,
                          "S.cerevisiae_tRNA-PHE M10740/1-73",
                          "(((((((..((((........)))).((((.........)))).....(((((.......))))))))))))."_wuss51);
        fout.emplace_back("UUGGAGUACACAACCUGUACACUCUUUC"_rna4,
                          "example",
                          "..(((((..(((...)))..)))))..."_wuss51);
    }

    seqan3::structure_file_input sf{tmp_file}; // Vienna with RNA sequences assumed, use regular std::ifstream as stream

    seqan3::structure_file_input fin{std::istringstream{input}, seqan3::format_vienna{}};
    //                          ^ no need to specify the template arguments

    std::filesystem::remove(tmp_file);
}
