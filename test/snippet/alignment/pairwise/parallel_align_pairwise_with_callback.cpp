// SPDX-FileCopyrightText: 2006-2025 Knut Reinert & Freie Universität Berlin
// SPDX-FileCopyrightText: 2016-2025 Knut Reinert & MPI für molekulare Genetik
// SPDX-License-Identifier: CC0-1.0

#include <mutex>
#include <vector>

#include <seqan3/alignment/configuration/all.hpp>
#include <seqan3/alignment/pairwise/align_pairwise.hpp>
#include <seqan3/alphabet/nucleotide/dna4.hpp>
#include <seqan3/core/debug_stream.hpp>

int main()
{
    // Generate some sequences.
    using namespace seqan3::literals;
    using sequence_pair_t = std::pair<seqan3::dna4_vector, seqan3::dna4_vector>;
    std::vector<sequence_pair_t> sequences{100, {"AGTGCTACG"_dna4, "ACGTGCGACTAG"_dna4}};

    // Use edit distance with 4 threads.
    auto const alignment_config =
        seqan3::align_cfg::method_global{} | seqan3::align_cfg::edit_scheme | seqan3::align_cfg::parallel{4};

    // Compute the alignments in parallel and output them in order based on the input.
    for (auto && result : seqan3::align_pairwise(sequences, alignment_config))
        seqan3::debug_stream << result << '\n';

    seqan3::debug_stream << '\n';
    // prints:
    // [id: 0 score: -4]
    // [id: 1 score: -4]
    // [id: 2 score: -4]
    // [id: 3 score: -4]
    // [id: 4 score: -4]
    // [id: 5 score: -4]
    // ...
    // [id: 98 score: -4]
    // [id: 99 score: -4]

    // Compute the alignments in parallel and output them unordered using the callback (order is not deterministic).
    std::mutex write_to_debug_stream{}; // Need mutex to synchronise the output.
    auto const alignment_config_with_callback =
        alignment_config
        | seqan3::align_cfg::on_result{[&](auto && result)
                                       {
                                           std::lock_guard sync{write_to_debug_stream}; // critical section
                                           seqan3::debug_stream << result << '\n';
                                       }};
    seqan3::align_pairwise(sequences, alignment_config_with_callback); // seqan3::align_pairwise is now declared void.

    // might print:
    // [id: 0 score: -4]
    // [id: 1 score: -4]
    // [id: 2 score: -4]
    // [id: 6 score: -4]
    // [id: 7 score: -4]
    // [id: 3 score: -4]
    // ...
    // [id: 99 score: -4]
    // [id: 92 score: -4]
}
