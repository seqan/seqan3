#include <string>
#include <vector>

#include <seqan3/alphabet/quality/phred42.hpp>
#include <seqan3/alphabet/views/to_char.hpp>
#include <seqan3/alphabet/views/trim_quality.hpp>

using namespace seqan3::literals;

int main()
{
    std::vector<seqan3::phred42> vec{"II?5+"_phred42};

    // trim by phred_value
    auto v1 = vec | seqan3::views::trim_quality(20u);                            // == ['I','I','?','5']

    // trim by quality character
    auto v2 = vec | seqan3::views::trim_quality('I'_phred42);                    // == ['I','I']

    // function syntax
    auto v3 = seqan3::views::trim_quality(vec, '5'_phred42);                     // == ['I','I','?','5']

    // combinability
    auto v4 = seqan3::views::trim_quality(vec, 20u) | seqan3::views::to_char;    // == "II?5"
}
