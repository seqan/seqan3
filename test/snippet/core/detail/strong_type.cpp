#include <vector>

#include <seqan3/core/detail/strong_type.hpp>
#include <seqan3/std/concept/range.hpp>
#include <seqan3/alphabet/nucleotide/dna4.hpp>

using namespace seqan3;

using namespace seqan3::literal;

std::vector<dna4> my_range = "ACGTT"_dna4;

namespace seqan3::detail
{
    template <forward_range_concept fwd_rng_type>
        bool do_find(fwd_rng_type &, int const, int const) { return true; }
}  // namespace seqan3::detail

//! [usage]
template <forward_range_concept fwd_rng_type>
bool search(fwd_rng_type & rng, unsigned const w, unsigned const e)
{
    // do something
    return detail::do_find(rng, w, e);
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
