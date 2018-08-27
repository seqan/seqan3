#include <seqan3/range/view/convert.hpp>
#include <seqan3/alphabet/nucleotide/dna15.hpp>
#include <seqan3/alphabet/nucleotide/dna5.hpp>
#include <seqan3/std/view/reverse.hpp>

using namespace seqan3;
using namespace seqan3::literal;

int main()
{

{
//! [int_to_bool]
// convert from int to bool
std::vector<int>  vec{7, 5, 0, 5, 0, 0, 4, 8, -3};

// pipe notation
auto v = vec | view::convert<bool>; // == [1, 1, 0, 1, 0, 0, 1, 1, 1];

// function notation and immediate conversion to vector again
std::vector<bool> v2(view::convert<bool>(vec));

// combinability
auto v3 = vec | view::convert<bool> | view::reverse; // == [1, 1, 1, 0, 0, 1, 0, 1, 1];
//! [int_to_bool]
(void) v;
(void) v2;
(void) v3;
}

{
//! [15_to_5]
dna15_vector vec2{"ACYGTN"_dna15};
auto v4 = vec2 | view::convert<dna5>; // == "ACNGTN"_dna5
//! [15_to_5]
(void) v4;
}

}
