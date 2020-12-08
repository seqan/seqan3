#include <seqan3/std/algorithm>

#include <seqan3/alphabet/quality/phred94.hpp>
#include <seqan3/core/debug_stream.hpp>

int main()
{
    using seqan3::operator""_phred94;

    // directly assign to a std::vector<phred94> using a string literal
    std::vector<seqan3::phred94> qual_vec = "###!"_phred94;

    // This is the same as a sequence of char literals:
    std::vector<seqan3::phred94> qual_vec2 = {'#'_phred94, '#'_phred94, '#'_phred94, '!'_phred94};

    seqan3::debug_stream << std::ranges::equal(qual_vec, qual_vec2) << '\n'; // prints 1 (true)
}
