#include <seqan3/core/platform.hpp>

#if !SEQAN3_WORKAROUND_GCC_96070
//![snippet]
#include <seqan3/std/iterator>
#include <seqan3/std/ranges>
#include <sstream>

#include <seqan3/io/sequence_file/input.hpp>
#include <seqan3/io/sequence_file/output.hpp>

auto input = R"(@TEST1
ACGT
+
##!#
@Test2
AGGCTGA
+
##!#!!!
@Test3
GGAGTATAATATATATATATATAT
+
##!###!###!###!###!###!#)";

int main()
{
    // minimum_average_quality_filter and minimum_sequence_length_filter need to be implemented first
    auto minimum_sequence_length_filter = std::views::filter([] (auto rec)
    {
        return std::ranges::distance(rec.sequence()) >= 50;
    });

    auto minimum_average_quality_filter = std::views::filter([] (auto const & record)
    {
        double qual_sum{0}; // summation of the qualities
        for (auto chr : record.base_qualities())
            qual_sum += seqan3::to_phred(chr);

        // check if average quality is greater than 20.
        return qual_sum / (std::ranges::distance(record.base_qualities())) >= 20;
    });

    auto input_file = seqan3::sequence_file_input{std::istringstream{input}, seqan3::format_fastq{}};
    input_file | minimum_average_quality_filter
               | minimum_sequence_length_filter
               | std::views::take(3)
               | seqan3::sequence_file_output{std::ostringstream{}, seqan3::format_fasta{}};
}
//![snippet]
#endif // !SEQAN3_WORKAROUND_GCC_96070
