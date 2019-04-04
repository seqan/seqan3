// -----------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2019, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2019, Knut Reinert & MPI für molekulare Genetik
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at: https://github.com/seqan/seqan3/blob/master/LICENSE
// -----------------------------------------------------------------------------------------------------

// written for the purpose of the Bch Softwarepraktikum, early 2019 FU Berlin
// @author: Clemens Cords <clemenscords@fu-berlin.de>

#include <cstring>
#include <random>

#include <benchmark/benchmark.h>

#include <seqan3/io/stream/debug_stream.hpp>
#include <seqan3/search/fm_index/all.hpp>
#include <seqan3/alphabet/all.hpp>

using namespace seqan3;

const std::size_t LENGTH = 50;
const unsigned int SEED = 1234;

// ============================================================================
// generate sequence deterministically based on fixed const seed
// ============================================================================
template<Alphabet alphabet_t>
std::vector<alphabet_t> generate_test_sequence(const std::size_t length)
{
    const size_t alphabet_size = alphabet_t::value_size;

    std::mt19937 rng_engine;
    rng_engine.seed(SEED);
    std::uniform_int_distribution distribution{0, static_cast<int>(alphabet_size-1)};

    std::vector<alphabet_t> text{};

    for (std::size_t i = 0; i < length; ++i)
    {
        alphabet_t s{};
        s.assign_rank(distribution(rng_engine));
        text.push_back(s);
    }

    return text;
}

// ============================================================================
// "overload" for chars
// ============================================================================
std::vector<char> generate_test_sequence_char(const std::size_t length)
{
    std::mt19937 rng_engine;
    rng_engine.seed(SEED);
    std::uniform_int_distribution distribution{0, 255};

    std::vector<char> text{};

    for (std::size_t i = 0; i < length; ++i)
    {
        char c = static_cast<char>(distribution(rng_engine));
        text.push_back(c);
    }

    return text;
}

// ============================================================================
// construct as many fm indexes as possible
// ============================================================================
template<Alphabet alphabet_t>
void construct_fm_index(benchmark::State & state)
{
    const auto sequence = generate_test_sequence<alphabet_t>(LENGTH);
    fm_index<std::vector<alphabet_t>> index{};

    for (auto _ : state)
    {
        index = fm_index{sequence};
    }
}

// ============================================================================
// construct as many bi fm indexes as possible
// ============================================================================
template<Alphabet alphabet_t>
void construct_bi_fm_index(benchmark::State & state)
{
    const auto sequence = generate_test_sequence<alphabet_t>(LENGTH);
    bi_fm_index<std::vector<alphabet_t>> index{};

    for (auto _ : state)
    {
        index = bi_fm_index{sequence};
    }
}

// ============================================================================
// comparison with 256 char alphabet
// ============================================================================
void construct_fm_index_char(benchmark::State & state)
{
    // hardcode because char does not model alphabet concept
    const auto sequence = generate_test_sequence_char(LENGTH);
    fm_index<std::vector<char>> index{};

    for (auto _ : state)
    {
        index = fm_index{sequence};
    }
}

void construct_bi_fm_index_char(benchmark::State & state)
{
    const auto sequence = generate_test_sequence_char(LENGTH);
    bi_fm_index<std::vector<char>> index{};

    for (auto _ : state)
    {
        index = bi_fm_index{sequence};
    }
}

BENCHMARK_TEMPLATE(construct_fm_index, dna4);
BENCHMARK_TEMPLATE(construct_fm_index, aa27);
BENCHMARK(construct_fm_index_char);
BENCHMARK_TEMPLATE(construct_bi_fm_index, dna4);
BENCHMARK_TEMPLATE(construct_bi_fm_index, aa27);
BENCHMARK(construct_bi_fm_index_char);

BENCHMARK_MAIN();