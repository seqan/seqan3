#include <seqan3/alphabet/quality/phred68legacy.hpp>
#include <seqan3/core/debug_stream.hpp>
#include <seqan3/std/algorithm>

int main()
{
    using seqan3::operator""_phred68legacy;

    // directly assign to a std::vector<phred68legacy> using a string literal
    std::vector<seqan3::phred68legacy> qual_vec = "###!"_phred68legacy;

    // This is the same as a sequence of char literals:
    std::vector<seqan3::phred68legacy> qual_vec2 = {'#'_phred68legacy, '#'_phred68legacy,
                                                    '#'_phred68legacy, '!'_phred68legacy};

    seqan3::debug_stream << std::ranges::equal(qual_vec, qual_vec2) << '\n'; // prints 1 (true)
}
