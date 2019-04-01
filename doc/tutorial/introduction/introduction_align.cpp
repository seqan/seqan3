#include <seqan3/io/sequence_file/output.hpp>

//! [sequence_input_include]
#include <seqan3/io/sequence_file/input.hpp> // for sequence_file_input
#include <seqan3/io/stream/debug_stream.hpp> // for debug_stream
//! [sequence_input_include]

//! [alignment_include]
#include <tuple>                                                          // for std::make_pair
#include <seqan3/alignment/aligned_sequence/aligned_sequence_concept.hpp> // for alignment stream operator
#include <seqan3/alignment/pairwise/align_pairwise.hpp>                   // for align_pairwise
//! [alignment_include]

using namespace seqan3;

int main()
{
    auto tmp_dir = std::filesystem::temp_directory_path();
    std::string filename{tmp_dir/"seq.fasta"};
    {
        // Create a /tmp/my.fasta file.
        sequence_file_output file_out{filename};
        file_out.emplace_back("ACGTGATG"_dna4, std::string{"seq1"});
        file_out.emplace_back("AGTGATACT"_dna4, std::string{"seq2"});
    }

//! [sequence_input]
    // Initialise a file input object with a FastA file.
    sequence_file_input file_in{filename}; // filename: "seq.fasta"

    // Retrieve the sequences and ids.
    for (auto &[seq, id, qual] : file_in)
    {
        debug_stream << "ID:  " << id << '\n';
        debug_stream << "SEQ: " << seq << '\n';
        debug_stream << "EMPTY QUAL." << qual << '\n'; // qual is empty for FastA files
    }
//! [sequence_input]

    std::vector<seqan3::dna5_vector> sequences{"ACGTGATG"_dna5, "AGTGATACT"_dna5};
//! [alignment]
    // Call a pairwise alignment with edit distance and traceback.
    for (auto && res : align_pairwise(std::tie(sequences[0], sequences[1]),
                                      align_cfg::edit | align_cfg::result{with_alignment}))
    {
        // Print the resulting score and the alignment.
        debug_stream << res.score() << '\n';              // => -4
        debug_stream << res.alignment() << '\n';          // =>       0     .    :
                                                          //            ACGTGATG--
                                                          //            | |||||
                                                          //            A-GTGATACT
    }
//! [alignment]
    return 0;
}
