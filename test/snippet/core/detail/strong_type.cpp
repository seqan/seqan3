#include <vector>

#include <seqan3/core/detail/strong_type.hpp>
#include <seqan3/std/concept/range.hpp>
#include <seqan3/alphabet/nucleotide/dna4.hpp>

using namespace seqan3;
using namespace seqan3::literal;

std::vector<dna4> my_range = "ACGTT"_dna4;
//! [usage]
template <forward_range_concept fwd_rng_type>
bool search(fwd_rng_type & SEQAN3_DOXYGEN_ONLY(rng), unsigned const SEQAN3_DOXYGEN_ONLY(w), unsigned const SEQAN3_DOXYGEN_ONLY(e))
{
    // do something
    return true;
}

int main()
{
    // do something
    search(my_range, 4, 2);
    // do something
}
//! [usage]

void test()
{
//! [adding_skills]
struct error : detail::strong_type<unsigned, error, detail::strong_type_skill::decrement |
                                                    detail::strong_type_skill::increment>
{
    using detail::strong_type<unsigned, error, detail::strong_type_skill::decrement |
                                               detail::strong_type_skill::increment>::strong_type;
};
error e{4};
--e;
++e;
}
//! [adding_skills]
