#include <gtest/gtest.h>

#include <seqan3/range/view/trim.hpp>
#include <seqan3/range/view/to_char.hpp>
#include <seqan3/alphabet/quality/phred42.hpp>
#include <seqan3/alphabet/quality/aliases.hpp>

using namespace seqan3;

int main()
{

{
//! [phred42]
std::vector<phred42> vec{phred42{40}, phred42{40}, phred42{30}, phred42{20}, phred42{10}};

// trim by phred_value
auto v1 = vec | view::trim(20u);                        // == ['I','I','?','5']

// trim by quality character
auto v2 = vec | view::trim(phred42{40});             // == ['I','I']

// function syntax
auto v3 = view::trim(vec, 20u);                         // == ['I','I','?','5']

// combinability
std::string v4 = view::trim(vec, 20u) | view::to_char;  // == "II?5"
//! [phred42]
(void) v1;
(void) v2;
(void) v3;
}

{
//! [dna5q]
std::vector<dna5q> vec{{'A'_dna5, phred42{40}}, {'G'_dna5, phred42{40}}, {'G'_dna5, phred42{30}},
                       {'A'_dna5, phred42{20}}, {'T'_dna5, phred42{10}}};
std::vector<dna5q> cmp{{'A'_dna5, phred42{40}}, {'G'_dna5, phred42{40}}, {'G'_dna5, phred42{30}},
                       {'A'_dna5, phred42{20}}};

// trim by phred_value
auto v1 = vec | view::trim(20u);
assert(std::vector<dna5q>(v1) == cmp);

// trim by quality character; in this case the nucleotide part of the character is irrelevant
auto v2 = vec | view::trim(dna5q{'C'_dna5, phred42{20}});
assert(std::vector<dna5q>(v2) == cmp);

// combinability
std::string v3 = view::trim(vec, 20u) | view::to_char;
EXPECT_EQ("AGGA", v3);
//! [dna5q]
(void) v1;
(void) v2;
}

}
