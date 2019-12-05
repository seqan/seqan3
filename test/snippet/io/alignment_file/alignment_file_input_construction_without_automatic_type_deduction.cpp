#include <sstream>

#include <seqan3/io/alignment_file/input.hpp>
#include <seqan3/core/type_list/type_list.hpp>

auto input = R"(@HD	VN:1.6	SO:coordinate
r001	99	ref	7	30	8M2I4M1D3M	=	37	39	TTAGATAAAGGATACTG	*)";

int main()
{
    // The default types; you can adjust this list if you don't want to read all this data.
    using default_fields = seqan3::fields<seqan3::field::SEQ,
                                          seqan3::field::ID,
                                          seqan3::field::OFFSET,
                                          seqan3::field::REF_SEQ,
                                          seqan3::field::REF_ID,
                                          seqan3::field::REF_OFFSET,
                                          seqan3::field::ALIGNMENT,
                                          seqan3::field::MAPQ,
                                          seqan3::field::QUAL,
                                          seqan3::field::FLAG,
                                          seqan3::field::MATE,
                                          seqan3::field::TAGS,
                                          seqan3::field::EVALUE,
                                          seqan3::field::BIT_SCORE,
                                          seqan3::field::HEADER_PTR>;

                                                                // The expected format:
    using alignment_file_input_t = seqan3::alignment_file_input<seqan3::alignment_file_input_default_traits<>,
                                                                default_fields,
                                                                // Which formats are allowed:
                                                                seqan3::type_list<seqan3::format_sam>>;

    alignment_file_input_t fin{std::istringstream{input}, seqan3::format_sam{}};
}
