// SPDX-FileCopyrightText: 2006-2025 Knut Reinert & Freie Universität Berlin
// SPDX-FileCopyrightText: 2016-2025 Knut Reinert & MPI für molekulare Genetik
// SPDX-License-Identifier: CC0-1.0

#include <mutex>
#include <vector>

#include <seqan3/alignment/configuration/align_config_edit.hpp>
#include <seqan3/alignment/configuration/align_config_method.hpp>
#include <seqan3/alignment/configuration/align_config_on_result.hpp>
#include <seqan3/alignment/configuration/align_config_parallel.hpp>
#include <seqan3/alignment/pairwise/align_pairwise.hpp>
#include <seqan3/alphabet/nucleotide/dna4.hpp>
#include <seqan3/core/debug_stream.hpp>

int main()
{
    // Generate some sequences.
    using seqan3::operator""_dna4;
    using sequence_pair_t = std::pair<seqan3::dna4_vector, seqan3::dna4_vector>;
    std::vector<sequence_pair_t> sequences{100, {"AGTGCTACG"_dna4, "ACGTGCGACTAG"_dna4}};

    std::mutex write_to_debug_stream{}; // Need mutex to synchronise the output.

    // Use edit distance with 4 threads.
    auto const alignment_config =
        seqan3::align_cfg::method_global{} | seqan3::align_cfg::edit_scheme | seqan3::align_cfg::parallel{4}
        | seqan3::align_cfg::on_result{[&](auto && result)
                                       {
                                           std::lock_guard sync{write_to_debug_stream}; // critical section
                                           seqan3::debug_stream << result << '\n';
                                       }};

    // Compute the alignments in parallel, and output them unordered using the callback (order is not deterministic).
    seqan3::align_pairwise(sequences, alignment_config); // seqan3::align_pairwise is now declared void.
}
