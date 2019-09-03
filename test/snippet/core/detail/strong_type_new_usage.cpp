#include <vector>

#include <seqan3/alphabet/nucleotide/dna4.hpp>
#include <seqan3/core/detail/strong_type.hpp>
#include <seqan3/std/ranges>

struct error : seqan3::detail::strong_type<unsigned, error>
{
    using seqan3::detail::strong_type<unsigned, error>::strong_type;
};

struct window_size : seqan3::detail::strong_type<unsigned, window_size>
{
    using seqan3::detail::strong_type<unsigned, window_size>::strong_type;
};

namespace seqan3::detail
{
template <std::ranges::forward_range fwd_rng_type>
    bool do_find(fwd_rng_type &, int const, int const) { return true; }
}  // namespace seqan3::detail

template <std::ranges::forward_range fwd_rng_type>
    bool search(fwd_rng_type & rng, window_size const w, error const e)
{
    // do something
    return seqan3::detail::do_find(rng, w.get(), e.get());
}

int main()
{
    using seqan3::operator""_dna4;

    std::vector<seqan3::dna4> my_range = "ACGTT"_dna4;
    // do something
    search(my_range, window_size{4}, error{2});
    // do something
    return 0;
}
