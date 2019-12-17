#include <sstream>

#include <seqan3/core/debug_stream.hpp>
#include <seqan3/core/type_list/type_list.hpp>
#include <seqan3/io/alignment_file/input.hpp>
#include <seqan3/std/ranges>

auto sam_file_raw = R"(@HD	VN:1.6	SO:coordinate	GO:none
@SQ	SN:ref	LN:45
r001	99	ref	7	30	8M2I4M1D3M	=	37	39	TTAGATAAAGGATACTG	*
r003	0	ref	29	30	5S6M	*	0	0	GCCTAAGCTAA	*	SA:Z:ref,29,-,6H5M,17,0;
r003	2064	ref	29	17	6H5M	*	0	0	TAGGC	*	SA:Z:ref,9,+,5S6M,30,1;
r001	147	ref	237	30	9M	=	7	-39	CAGCGGCAT	*	NM:i:1
)";

int main()
{
    using seqan3::get;

    seqan3::alignment_file_input fin{std::istringstream{sam_file_raw}, seqan3::format_sam{}};

    auto minimum_length10_filter = std::views::filter([] (auto const & rec)
    {
        return std::ranges::size(get<seqan3::field::seq>(rec)) >= 10;
    });

    for (auto & rec : fin | minimum_length10_filter) // only records with sequence length >= 10 will "appear"
        seqan3::debug_stream << get<seqan3::field::id>(rec) << '\n';
}
