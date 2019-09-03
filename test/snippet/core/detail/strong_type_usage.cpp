#include <vector>

#include <seqan3/alphabet/nucleotide/dna4.hpp>
#include <seqan3/std/ranges>

namespace seqan3::detail
{
template <std::ranges::forward_range fwd_rng_type>
bool do_find(fwd_rng_type &, int const, int const) { return true; }
}  // namespace seqan3::detail

template <std::ranges::forward_range fwd_rng_type>
bool search(fwd_rng_type & rng, unsigned const w, unsigned const e)
{
    // do something
    return seqan3::detail::do_find(rng, w, e);
}

int main()
{
    using seqan3::operator""_dna4;

    std::vector<seqan3::dna4> my_range = "ACGTT"_dna4;
    // do something
    search(my_range, 4, 2);
    // do something
    return 0;
}
