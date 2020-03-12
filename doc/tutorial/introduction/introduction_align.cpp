//! [alignment]
#include <vector>                                         // for std::vector
#include <tuple>                                          // for std::tie

#include <seqan3/alignment/aligned_sequence/all.hpp>      // for alignment stream operator <<
#include <seqan3/alignment/configuration/all.hpp>         // for all configs in the seqan3::align_cfg namespace
#include <seqan3/alignment/pairwise/all.hpp>              // for seqan3::align_pairwise
#include <seqan3/alphabet/nucleotide/dna5.hpp>            // for dna5 datastrucutre
#include <seqan3/argument_parser/all.hpp>                 // for argument_parser
#include <seqan3/core/debug_stream.hpp>                   // for debug_stream
#include <seqan3/io/sequence_file/all.hpp>                // for sequence_file_input and sequence_file_output
#include <seqan3/std/filesystem>                          // for tmp_dir

int main()
{
    auto tmp_dir = std::filesystem::temp_directory_path();
    std::string filename{tmp_dir/"seq.fasta"};
    {
        // Create a /tmp/my.fasta file.
        seqan3::sequence_file_output file_out{filename};

        using seqan3::operator""_dna5;

        file_out.emplace_back("ACGTGATG"_dna5, std::string{"seq1"});
        file_out.emplace_back("AGTGATACT"_dna5, std::string{"seq2"});
    }

    // Initialise a file input object and a vector
    seqan3::sequence_file_input file_in{filename};
    std::vector<seqan3::dna5_vector> sequences;

    for (auto & [ seq, id, qual ] : file_in)
    {
        sequences.push_back(seq);
    }

    // Call a pairwise alignment with edit distance and traceback.
    for (auto && res : align_pairwise(std::tie(sequences[0], sequences[1]),
                                      seqan3::align_cfg::edit | seqan3::align_cfg::result{seqan3::with_alignment}))
    {
        // Print the resulting score and the alignment.
        seqan3::debug_stream << res.score() << '\n';      // => -4
        seqan3::debug_stream << res.alignment() << '\n';  // =>       0     .    :
                                                          //            ACGTGATG--
                                                          //            | |||||
                                                          //            A-GTGATACT
    }
    std::filesystem::remove(filename);
    return 0;
}
//! [alignment]
