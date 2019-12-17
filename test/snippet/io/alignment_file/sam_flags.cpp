#include <iostream>
#include <sstream>

#include <seqan3/io/alignment_file/all.hpp>

auto sam_file_raw = R"(@HD	VN:1.6	SO:coordinate	GO:none
@SQ	SN:ref	LN:45
r001	99	ref	7	30	8M2I4M1D3M	=	37	39	TTAGATAAAGGATACTG	!!!!!!!!!!!!!!!!!
r003	0	ref	29	30	5S6M	*	0	0	GCCTAAGCTAA	!!!!!!!!!!!	SA:Z:ref,29,-,6H5M,17,0;
r003	4	*	29	17	*	*	0	0	TAGGC	@@@@@	SA:Z:ref,9,+,5S6M,30,1;
r001	147	ref	237	30	9M	=	7	-39	CAGCGGCAT	!!!!!!!!!	NM:i:1
)";

int main()
{
    seqan3::alignment_file_input fin{std::istringstream{sam_file_raw}, seqan3::format_sam{}};

    for (auto & rec : fin)
    {
        // Check if a certain flag value (bit) is set:
        if (static_cast<bool>(seqan3::get<seqan3::field::flag>(rec) & seqan3::sam_flag::unmapped))
            std::cout << "Read " << seqan3::get<seqan3::field::id>(rec) << " is unmapped\n";

        if (seqan3::get<seqan3::field::qual>(rec)[0] < seqan3::assign_char_to('@', seqan3::phred42{})) // low quality
        {
            // Set a flag value (bit):
            seqan3::get<seqan3::field::flag>(rec) &= seqan3::sam_flag::failed_filter;
            // Note that this does not affect other flag values (bits),
            // e.g. `seqan3::get<seqan3::field::flag>(rec) & seqan3::sam_flag::unmapped` may still be true
        }

        // Unset a flag value (bit):
        seqan3::get<seqan3::field::flag>(rec) &= ~seqan3::sam_flag::duplicate; // not marked as a duplicate anymore
    }
}
