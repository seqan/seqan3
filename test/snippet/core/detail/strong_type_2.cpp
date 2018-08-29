#include <vector>

#include <seqan3/core/detail/strong_type.hpp>
#include <seqan3/std/ranges>
#include <seqan3/alphabet/nucleotide/dna4.hpp>

using namespace seqan3;

//! [error_window]
struct error : detail::strong_type<unsigned, error>
{
    using detail::strong_type<unsigned, error>::strong_type;
};

struct window_size : detail::strong_type<unsigned, window_size>
{
    using detail::strong_type<unsigned, window_size>::strong_type;
};
//! [error_window]

using namespace seqan3::literal;

namespace seqan3::detail
{
template <std::ranges::ForwardRange fwd_rng_type>
    bool do_find(fwd_rng_type &, int const, int const) { return true; }
}  // namespace seqan3::detail

std::vector<dna4> my_range = "ACGTT"_dna4;
//! [new_usage]
template <std::ranges::ForwardRange fwd_rng_type>
    bool search(fwd_rng_type & rng, window_size const w, error const e)
{
    // do something
    return detail::do_find(rng, w.get(), e.get());
}

int main()
{
    // do something
    search(my_range, window_size{4}, error{2});
    // do something
    return 0;
}
//! [new_usage]
