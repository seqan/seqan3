#include <seqan3/range/view/deep.hpp>
#include <seqan3/range/view/take.hpp>
#include <seqan3/std/view/reverse.hpp>
#include <seqan3/alphabet/nucleotide/dna5.hpp>

namespace seqan3::view
{
inline auto const deep_reverse = deep{view::reverse};
inline auto const deep_take = deep{view::take};
inline auto const deep_take1 = deep{view::take(1)};
}

using namespace seqan3;
using namespace seqan3::literal;

int main()
{
{
std::vector<dna5_vector> foo{"AAATTT"_dna5, "CCCGGG"_dna5};

auto r = foo | view::reverse;             // == [ [C,C,C,G,G,G], [A,A,A,T,T,T] ]

auto d = foo | view::deep{view::reverse}; // == [ [T,T,T,A,A,A], [G,G,G,C,C,C] ]

auto e = foo | view::deep_reverse;                // == [ [T,T,T,A,A,A], [G,G,G,C,C,C] ]

(void) r;
(void) d;
(void) e;
#if 0 // Create a copy of the code in a comment to include the namespace declaration.
//! [no_param]
std::vector<dna5_vector> foo{"AAATTT"_dna5, "CCCGGG"_dna5};

auto r = foo | view::reverse;             // == [ [C,C,C,G,G,G], [A,A,A,T,T,T] ]

auto d = foo | view::deep{view::reverse}; // == [ [T,T,T,A,A,A], [G,G,G,C,C,C] ]

// You can also create a permanent alias:
namespace view
{
inline auto const deep_reverse = deep{view::reverse};
}

auto e = foo | view::deep_reverse;        // == [ [T,T,T,A,A,A], [G,G,G,C,C,C] ]
//! [no_param]
#endif
}

{
std::vector<dna5_vector> foo{"AAATTT"_dna5, "CCCGGG"_dna5};

auto t = foo | view::take(1);             // == [ [A,A,A,T,T,T] ]

auto d = foo | view::deep{view::take}(1); // == [ [A], [C] ]
// constructor arguments passed via {} and arguments to underlying view passed via ()

auto e = foo | view::deep_take(1);                // == [ [A], [C] ]

(void) t;
(void) d;
(void) e;
#if 0 // Copied code for documentation.
//! [with_param]
std::vector<dna5_vector> foo{"AAATTT"_dna5, "CCCGGG"_dna5};

auto t = foo | ranges::view::take(1);             // == [ [A,A,A,T,T,T] ]

auto d = foo | view::deep{ranges::view::take}(1); // == [ [A], [C] ]
// constructor arguments passed via {} and arguments to underlying view passed via ()

// In this case especially, an alias improves readability:
namespace view
{
inline auto const deep_take = deep{ranges::view::take};
}

auto e = foo | view::deep_take(1);                // == [ [A], [C] ]
//! [with_param]
#endif

//! [pass_ref]
int i = 7;
auto f = foo | view::deep_take(i);
//! [pass_ref]
(void) f;
}

{
std::vector<dna5_vector> foo{"AAATTT"_dna5, "CCCGGG"_dna5};

auto t = foo | view::take(1);             // == [ [A,A,A,T,T,T] ]

auto d = foo | view::deep{view::take(1)}; // == [ [A], [C] ]
// constructor arguments passed via {} and arguments to underlying view hardcoded inside

auto e = foo | view::deep_take1;                  // == [ [A], [C] ]

(void) t;
(void) d;
(void) e;
#if 0 //Copied code for documentation.
//! [wrap_args]
std::vector<dna5_vector> foo{"AAATTT"_dna5, "CCCGGG"_dna5};

auto t = foo | ranges::view::take(1);             // == [ [A,A,A,T,T,T] ]

auto d = foo | view::deep{ranges::view::take(1)}; // == [ [A], [C] ]
// constructor arguments passed via {} and arguments to underlying view hardcoded inside

// Or with an alias
namespace view
{
inline auto const deep_take1 = deep{ranges::view::take(1)};
}

auto e = foo | view::deep_take1;                  // == [ [A], [C] ]
//! [wrap_args]
#endif
}
}
