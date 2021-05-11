#include <seqan3/std/ranges>
#include <vector>

#include <seqan3/alphabet/nucleotide/dna4.hpp>

namespace detail
{

template <std::ranges::forward_range fwd_rng_type>
bool do_find(fwd_rng_type const &, uint8_t const, uint8_t const)
{
    return true;
}

}  // namespace detail

template <std::ranges::forward_range fwd_rng_type>
bool search(fwd_rng_type const & rng, uint8_t const window_size, uint8_t const error)
{
    return detail::do_find(rng, window_size, error);
}

int main()
{
    using namespace seqan3::literals;

    std::vector<seqan3::dna4> range = "ACGTT"_dna4;
    search(range, 4u, 2u);

    return 0;
}
