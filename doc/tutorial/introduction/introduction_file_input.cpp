#include <seqan3/std/filesystem>                // to create temp directory
#include <seqan3/alphabet/nucleotide/dna4.hpp>  // to create right datastructure in tmp file
#include <seqan3/io/sequence_file/output.hpp>   // to create tmp file

// Creating a tmp FastA file with sequences.
// It's necessary, because some review processes in GitHub don't allow using
// an alreadys existing file. So you can use the whole code from this file as
// an example or the collection of snippets, shown in the introduction.

struct write_file_dummy_struct
{
    std::filesystem::path const file_path = std::filesystem::temp_directory_path()/"seq.fasta";

    write_file_dummy_struct()
    {

auto file_raw = R"//![fasta_file](
>seq1
ACGTGATG
>seq2
AGTGATACT
)//![fasta_file]";

        std::ofstream file{file_path};
        std::string str{file_raw};
        file << str.substr(1); // skip first newline
    }

    ~write_file_dummy_struct()
    {
        std::error_code ec{};
        std::filesystem::remove(file_path, ec);

        if (ec)
            seqan3::debug_stream << "[WARNING] Could not delete " << file_path << ". " << ec.message() << '\n';
    }
};

write_file_dummy_struct go{};

//! [fileinput]
#include <string>
#include <seqan3/core/debug_stream.hpp>         // for debug_stream
#include <seqan3/io/sequence_file/input.hpp>    // for sequence_file_input

int main ()
{
    // Initialise a file input object with a FastA file.
    std::filesystem::path tmp_dir = std::filesystem::temp_directory_path(); // get the temp directory

    seqan3::sequence_file_input file_in{tmp_dir/"seq.fasta"};

    // Retrieve the sequences and ids.
    for (auto &[seq, id, qual] : file_in)
    {
        seqan3::debug_stream << "ID:     " << id << '\n';
        seqan3::debug_stream << "SEQ:    " << seq << '\n';
        seqan3::debug_stream << "Empty Qual." << qual << '\n';  // qual es empty for FastAfiles
    }

    return 0;
}
//! [fileinput]
