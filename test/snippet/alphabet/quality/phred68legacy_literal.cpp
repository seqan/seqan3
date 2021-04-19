#include <seqan3/alphabet/quality/phred68legacy.hpp>
#include <seqan3/core/debug_stream.hpp>
#include <seqan3/std/algorithm>

int main()
{
    using seqan3::operator""_phred68solexa;

    // directly assign to a std::vector<phred68solexa> using a string literal
    std::vector<seqan3::phred68solexa> qual_vec = "###!"_phred68solexa;

    // This is the same as a sequence of char literals:
    std::vector<seqan3::phred68solexa> qual_vec2 = {'#'_phred68solexa, '#'_phred68solexa,
                                                    '#'_phred68solexa, '!'_phred68solexa};

    seqan3::debug_stream << std::ranges::equal(qual_vec, qual_vec2) << '\n'; // prints 1 (true)
}
