// SPDX-FileCopyrightText: 2006-2025 Knut Reinert & Freie Universität Berlin
// SPDX-FileCopyrightText: 2016-2025 Knut Reinert & MPI für molekulare Genetik
// SPDX-License-Identifier: BSD-3-Clause

/* Copied and adjusted from https://raw.githubusercontent.com/kloetzl/libdna/master/bench2/revcomp.cxx
 * Credits go to Fabian Klötzl (@kloetzl - https://github.com/kloetzl)
 */

#include <benchmark/benchmark.h>

#include <random>

#include <seqan3/alphabet/nucleotide/dna4.hpp>
#include <seqan3/alphabet/views/char_to.hpp>
#include <seqan3/alphabet/views/complement.hpp>
#include <seqan3/test/performance/simd_dna4.hpp>

enum class tag
{
    revcomp_dna4_inline,
    seqan3_dna4,
    seqan3_dna4_vector,
    seqan3_dna4_simd,
    seqan3_dna4_simd_vector
};

class allocator
{
public:
    allocator() = delete;
    allocator(allocator const &) = delete;
    allocator & operator=(allocator const &) = delete;
    allocator(allocator &&) = delete;
    allocator & operator=(allocator &&) = delete;
    ~allocator()
    {
        std::free(forward);
        std::free(reverse);
    }

    allocator(size_t const length)
    {
        forward = static_cast<char *>(std::malloc(length + 1));
        reverse = static_cast<char *>(std::malloc(length + 1));
    }

    auto get() const
    {
        return std::make_tuple(forward, reverse);
    }

private:
    char * forward{nullptr};
    char * reverse{nullptr};
};

static void generate_random_dna4_char_string(char * dest, size_t const length)
{
    constexpr std::array<char, 8> dna4_chars{'A', 'C', 'G', 'T', 'a', 'c', 'g', 't'};
    std::mt19937_64 random_engine{1729u};
    std::uniform_int_distribution<size_t> random_index(0u, dna4_chars.size() - 1u);
    auto random_dna4_char = [&dna4_chars, &random_engine, &random_index]()
    {
        return dna4_chars[random_index(random_engine)];
    };

    std::ranges::generate_n(dest, length, random_dna4_char);
    dest[length - 1u] = '\0';
}

static constexpr char * revcomp_dna4_inline(char const * const begin, size_t const length, char * const dest)
{
    for (size_t i = 0; i < length; ++i)
    {
        char c = begin[length - 1 - i];

        dest[i] = c ^= c & 2 ? 4 : 21;
    }

    return dest + length;
}

static constexpr auto seqan3_dna4(std::string_view sv)
{
    return sv | seqan3::views::char_to<seqan3::dna4> | std::views::reverse | seqan3::views::complement;
}

static void seqan3_dna4_vector(std::string_view sv, std::vector<seqan3::dna4> & dest)
{
    auto revcomp = sv | seqan3::views::char_to<seqan3::dna4> | std::views::reverse | seqan3::views::complement;
    std::ranges::copy(revcomp, dest.begin());
}

static constexpr auto seqan3_dna4_simd(std::string_view sv)
{
    return sv | seqan3::views::char_to<seqan3::simd_dna4> | std::views::reverse | seqan3::views::complement;
}

static void seqan3_dna4_simd_vector(std::string_view sv, std::vector<seqan3::simd_dna4> & dest)
{
    auto revcomp = sv | seqan3::views::char_to<seqan3::simd_dna4> | std::views::reverse | seqan3::views::complement;
    std::ranges::copy(revcomp, dest.begin());
}

template <tag id>
void complement(benchmark::State & state)
{
    constexpr size_t length{1'000'003};
    allocator alloc{length};
    auto [forward, revcomp] = alloc.get();
    generate_random_dna4_char_string(forward, length);

    auto sv = std::string_view{forward};
    using alphabet_t = std::conditional_t<id == tag::seqan3_dna4_vector, seqan3::dna4, seqan3::simd_dna4>;
    std::vector<alphabet_t> vector(length);

    for (auto _ : state)
    {
        if constexpr (id == tag::revcomp_dna4_inline)
        {
            revcomp_dna4_inline(forward, length, revcomp);
            benchmark::DoNotOptimize(revcomp);
        }
        else if constexpr (id == tag::seqan3_dna4)
        {
            auto view = seqan3_dna4(sv);
            for (auto elem : view)
                benchmark::DoNotOptimize(elem);
        }
        else if constexpr (id == tag::seqan3_dna4_vector)
        {
            seqan3_dna4_vector(sv, vector);
            benchmark::DoNotOptimize(vector);
        }
        else if constexpr (id == tag::seqan3_dna4_simd)
        {
            auto view = seqan3_dna4_simd(sv);
            for (auto elem : view)
                benchmark::DoNotOptimize(elem);
        }
        else if constexpr (id == tag::seqan3_dna4_simd_vector)
        {
            seqan3_dna4_simd_vector(sv, vector);
            benchmark::DoNotOptimize(vector);
        }
        else
        {
            throw std::logic_error{"Invalid tag"};
        }
    }
}

BENCHMARK_TEMPLATE(complement, tag::revcomp_dna4_inline);
BENCHMARK_TEMPLATE(complement, tag::seqan3_dna4);
BENCHMARK_TEMPLATE(complement, tag::seqan3_dna4_vector);
BENCHMARK_TEMPLATE(complement, tag::seqan3_dna4_simd);
BENCHMARK_TEMPLATE(complement, tag::seqan3_dna4_simd_vector);
