#include <sstream>

#include <seqan3/io/sequence_file/input.hpp>
#include <seqan3/io/sequence_file/output.hpp>
#include <seqan3/range/views/persist.hpp>
#include <seqan3/std/iterator>
#include <seqan3/std/ranges>

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
    using seqan3::get;

    // minimum_average_quality_filter and minimum_sequence_length_filter need to be implemented first
    auto minimum_sequence_length_filter = std::views::filter([] (auto rec)
    {
        return std::ranges::distance(get<seqan3::field::seq>(rec)) >= 50;
    });

    auto minimum_average_quality_filter = std::views::filter([] (auto const & rec)
    {
        double qual_sum{0}; // summation of the qualities
        for (auto chr : get<seqan3::field::qual>(rec))
            qual_sum += chr.to_phred();

                           // check if average quality is greater than 20.
        return qual_sum / (std::ranges::distance(get<seqan3::field::qual>(rec))) >= 20;
    });

    seqan3::sequence_file_input{std::istringstream{input}, seqan3::format_fastq{}}
        | seqan3::views::persist
        | minimum_average_quality_filter
        | minimum_sequence_length_filter
        | std::views::take(3)
        | seqan3::sequence_file_output{std::ostringstream{}, seqan3::format_fasta{}};
}
