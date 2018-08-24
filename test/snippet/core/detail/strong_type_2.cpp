#include <vector>

#include <seqan3/core/detail/strong_type.hpp>
#include <seqan3/std/concept/range.hpp>
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
std::vector<dna4> my_range = "ACGTT"_dna4;
//! [new_usage]
template <forward_range_concept fwd_rng_type>
    bool search(fwd_rng_type & SEQAN3_DOXYGEN_ONLY(rng), window_size const SEQAN3_DOXYGEN_ONLY(w), error const SEQAN3_DOXYGEN_ONLY(e))
{
    // do something
    return true;
}

int main()
{
    // do something
    search(my_range, window_size{4}, error{2});
    // do something
}
//! [new_usage]
