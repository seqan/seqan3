#include <seqan3/core/platform.hpp>

#ifdef SEQAN3_DEPRECATED_310
#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wdeprecated-declarations"

#include <seqan3/std/ranges>
#include <string>

#include <seqan3/range/views/as_const.hpp>

template <std::ranges::random_access_range rng_t>
void foobar(rng_t const & range)
{
    range[0] = 'A';
}

int main()
{
    std::string s{"CCC"};
    auto v0 = std::views::all(s);
    foobar(v0); // this is valid and will update s to "ACC"
                // because const-ness of view does not protect elements

    auto v1 = seqan3::views::as_const(s);
    // foobar(v1); // this is invalid, views::as_const protects elements
}

#pragma GCC diagnostic pop
#endif // SEQAN3_DEPRECATED_310
