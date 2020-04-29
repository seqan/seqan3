#include <vector>                               // std::vector
#include <forward_list>                         // std::forward_list
#include <seqan3/core/type_traits/lazy.hpp>     // seqan3::lazy and seqan3::lazy_conditional_t
#include <seqan3/core/type_traits/range.hpp>    // seqan3::size_type_t
#include <seqan3/std/ranges>                    // std::ranges::input_range

template <std::ranges::input_range rng_t>
void foobar(rng_t && range)
{

#if 0
    // The following would fail to compile if rng_t is not sized,
    // because std::ranges::range_size_t<rngt_t> needs to be valid
    // (independent of whether the condition evaluates to true)
    using size_type = std::conditional_t<std::ranges::sized_range<rng_t>,
                                         std::ranges::range_size_t<rng_t>,
                                         void>;
#endif

    // This delays instantiation of std::ranges::range_size_t<rngt_t> until after the
    // conditional-decision is made:
    using size_type = seqan3::detail::lazy_conditional_t<std::ranges::sized_range<rng_t>,
                                                         seqan3::detail::lazy<seqan3::size_type_t, rng_t>,
                                                         void>;

    // DO SOMETHING with size_type
}

int main()
{
    foobar(std::vector<int>{});         // sized
    foobar(std::forward_list<int>{});   // not sized
}
