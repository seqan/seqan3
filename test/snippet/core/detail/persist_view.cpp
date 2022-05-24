#include <seqan3/alphabet/nucleotide/dna4.hpp>
#include <seqan3/alphabet/views/to_char.hpp>
#include <seqan3/core/detail/persist_view.hpp>

// P2415R2 makes our persist view superfluous.
// It's implemented in GCC 12 and can be detected by checking __cpp_lib_ranges.
// P2415R2 allows binding of rvalues with std::views::owning_view and also adapts
// std::views::all to return an owning_view when applicable.
#if defined(__cpp_lib_ranges) && (__cpp_lib_ranges < 202110L)
#    define SEQAN3_STL_HAS_OWNING_VIEW 0
#else
#    define SEQAN3_STL_HAS_OWNING_VIEW 1
#endif

int main()
{
    using namespace seqan3::literals;

    // explicitly create an l-value of our dna vector:
    auto vec = "ACGT"_dna4;
    [[maybe_unused]] auto v = vec | seqan3::views::to_char;

    // using seqan3::detail::persist you can bind the temporary directly:
#if !SEQAN3_STL_HAS_OWNING_VIEW
    [[maybe_unused]] auto v2 = "ACGT"_dna4 | seqan3::detail::persist | seqan3::views::to_char;
#else
    [[maybe_unused]] auto v2 = "ACGT"_dna4 | seqan3::views::to_char;
#endif

    // note that seqan3::detail::persist must follow immediately after the temporary,
    // thus the function notation might be more intuitive:
#if !SEQAN3_STL_HAS_OWNING_VIEW
    [[maybe_unused]] auto v3 = seqan3::detail::persist("ACGT"_dna4) | seqan3::views::to_char;
#else
    [[maybe_unused]] auto v3 = seqan3::views::to_char("ACGT"_dna4);
#endif
}
