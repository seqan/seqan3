// -----------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2019, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2019, Knut Reinert & MPI für molekulare Genetik
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at: https://github.com/seqan/seqan3/blob/master/LICENSE
// -----------------------------------------------------------------------------------------------------

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

//! [index_search_include]
#include <seqan3/search/algorithm/all.hpp> // for fm_index, search
//! [index_search_include]

using namespace seqan3;

int main()
{
    auto tmp_dir = std::filesystem::temp_directory_path();
    {
        // Create a /tmp/my.fasta file.
        sequence_file_output file_out{tmp_dir/"seq.fasta"};
        file_out.emplace_back("ACGTGATG"_dna4, std::string{"seq1"});
        file_out.emplace_back("AGTGATACT"_dna4, std::string{"seq2"});
    }

//! [sequence_input]
    // Initialise a file input object with a FastA file.
    sequence_file_input file_in{tmp_dir/"seq.fasta"};

    // Retrieve the sequences and ids.
    auto & seqs = get<field::SEQ>(file_in);
    auto & ids  = get<field::ID>(file_in);

    // Print the content of the retrieved vectors.
    debug_stream << seqs << '\n';   // => [ACGTGATG,AGTGATACT]
    debug_stream << ids << '\n';    // => [seq1,seq2]
//! [sequence_input]

//! [alignment]
    // Call a pairwise alignment with edit distance and traceback.
    for (auto && res : align_pairwise(std::make_pair(seqs[0], seqs[1]),
                                      align_cfg::edit | align_cfg::result{align_cfg::with_trace}))
    {
        // Print the resulting score and the alignment.
        debug_stream << res.get_score() << '\n';              // => -4
        debug_stream << res.get_alignment() << '\n';          // =>       0     .    :
                                                              //            ACGTGATG--
                                                              //            | |||||
                                                              //            A-GTGATACT
    }
//! [alignment]

//! [index_search]
    // Define a genomic sequence.
    std::vector<dna4> genome{"ATCGATCGAAGGCTAGCTAGCTAAGGGA"_dna4};

    // Build an FM index for the given sequence.
    fm_index index{genome};

    // Search the TAG motiv in the index and print the positions, where it is found.
    debug_stream << search(index, "TAG"_dna4) << '\n';   // => [13,17]
//! [index_search]
    std::filesystem::remove(tmp_dir/"seq.fasta");
}
