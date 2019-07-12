#include <seqan3/alignment/configuration/align_config_aligned_ends.hpp>
#include <seqan3/core/debug_stream.hpp>

int main()
{
    // Setup for overlap alignment.
    seqan3::align_cfg::aligned_ends overlap{seqan3::free_ends_all};

    // Setup for global alignment.
    seqan3::align_cfg::aligned_ends global{seqan3::free_ends_none};

    // Setup for semi-global alignment with free-end gaps in the first sequence.
    seqan3::align_cfg::aligned_ends semi_seq1{seqan3::free_ends_first};

    // Setup for semi-global alignment with free-end gaps in the second sequence.
    seqan3::align_cfg::aligned_ends semi_seq2{seqan3::free_ends_second};

    // Custom settings.
    seqan3::align_cfg::aligned_ends custom{seqan3::end_gaps{seqan3::front_end_first{std::true_type{}}, seqan3::front_end_second{std::true_type{}}}};
}
