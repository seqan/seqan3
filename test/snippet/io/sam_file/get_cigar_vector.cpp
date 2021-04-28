#include <seqan3/alphabet/gap/gapped.hpp>
#include <seqan3/alphabet/nucleotide/dna4.hpp>
#include <seqan3/core/debug_stream.hpp>
#include <seqan3/io/sam_file/detail/cigar.hpp>

int main()
{
    using aligned_t = std::vector<seqan3::gapped<seqan3::dna4>>;
    using namespace seqan3::literals;

    aligned_t ref{'A'_dna4, 'T'_dna4, 'G'_dna4, 'G'_dna4, seqan3::gap{}, seqan3::gap{},
                  'C'_dna4, 'G'_dna4, 'T'_dna4, 'A'_dna4, 'G'_dna4, 'A'_dna4, 'G'_dna4, 'C'_dna4};

    aligned_t query{'A'_dna4, 'T'_dna4, 'G'_dna4, 'C'_dna4, 'C'_dna4, 'C'_dna4, 'C'_dna4,
                    'G'_dna4, 'T'_dna4, 'T'_dna4, 'G'_dna4, seqan3::gap{}, seqan3::gap{}, 'C'_dna4};

    seqan3::debug_stream << seqan3::detail::get_cigar_vector(std::tie(ref, query));
}
