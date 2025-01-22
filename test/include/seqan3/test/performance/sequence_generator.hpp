// SPDX-FileCopyrightText: 2006-2025 Knut Reinert & Freie Universität Berlin
// SPDX-FileCopyrightText: 2016-2025 Knut Reinert & MPI für molekulare Genetik
// SPDX-License-Identifier: BSD-3-Clause

/*!\file
 * \brief Provides test utilities for seqan3::simd::simd_type types.
 * \author Marcel Ehrhardt <marcel.ehrhardt AT fu-berlin.de>
 * \author Lydia Buntrock <lydia.buntrock AT fu-berlin.de>
 */

#pragma once

#include <algorithm>
#include <random>
#include <ranges>

#include <seqan3/alphabet/concept.hpp>
#include <seqan3/test/seqan2.hpp>

#ifdef SEQAN3_HAS_SEQAN2
#    include <seqan/basic.h>
#    include <seqan/sequence.h>
#endif

namespace seqan3::test
{

template <typename sequence_t>
    requires (std::ranges::output_range<sequence_t, std::ranges::range_value_t<sequence_t>> &&
              requires (sequence_t & s) { {s.resize(1)}; })
class random_sequence_generator
{
public:
    //!\brief Stores the size of the random sequence.
    size_t size{};
    //!\brief Stores the variance of the random sequence.
    size_t size_variance{};

    random_sequence_generator() = default;                                              //!< Defaulted
    random_sequence_generator(random_sequence_generator const &) = default;             //!< Defaulted
    random_sequence_generator(random_sequence_generator &&) = default;                  //!< Defaulted
    random_sequence_generator & operator=(random_sequence_generator const &) = default; //!< Defaulted
    random_sequence_generator & operator=(random_sequence_generator &&) = default;      //!< Defaulted
    ~random_sequence_generator() = default;                                             //!< Defaulted

    /*!\brief Initialises a random sequence generator which generates sequences with a given mean size.
     * \param[in] size The size of the random sequence.
     * \param[in] size_variance The variance of the random sequence size (defaults to `0`).
     */
    explicit random_sequence_generator(size_t size, size_t size_variance = 0) : size{size}, size_variance{size_variance}
    {}

    /*!\name Element access
     * \{
     */

    /*!\brief Returns a random sequence for a set size (and variance).
     * \param[in,out] random_generator, e.g. std::mt19937 with a seed of 42.
     * \returns A generated random sequence.
     */
    template <typename generator_t>
    sequence_t operator()(generator_t && random_generator) const
    {
        using alphabet_t = std::ranges::range_value_t<sequence_t>;

        size_t max_alphabet_rank{};

        if constexpr (std::is_same_v<alphabet_t, uint64_t>)
            max_alphabet_rank = std::numeric_limits<size_t>::max();
        else
            max_alphabet_rank = seqan3::alphabet_size<alphabet_t> - 1ull;

        std::uniform_int_distribution<size_t> alphabet_rank_distribution(0ull, max_alphabet_rank);
        std::uniform_int_distribution<size_t> sequence_size_distribution(size - size_variance, size + size_variance);

        size_t const sequence_size = sequence_size_distribution(random_generator);
        sequence_t random_sequence;
        random_sequence.resize(sequence_size);
        std::ranges::generate(random_sequence,
                              [&]() -> alphabet_t
                              {
                                  // uint64_t is not an alphabet, see https://github.com/seqan/product_backlog/issues/200
                                  if constexpr (std::is_same_v<alphabet_t, uint64_t>)
                                      return alphabet_rank_distribution(random_generator);
                                  else
                                      return seqan3::assign_rank_to(alphabet_rank_distribution(random_generator),
                                                                    alphabet_t{});
                              });

        return random_sequence;
    }
    //!\}
};

template <typename alphabet_t>
auto generate_sequence(size_t const size = 500, size_t const size_variance = 0, size_t const seed = 0)
{
    seqan3::test::random_sequence_generator<std::vector<alphabet_t>> random_sequence_gen{size, size_variance};
    std::mt19937_64 random_engine{seed};

    return random_sequence_gen(random_engine);
}

template <arithmetic number_type>
auto generate_numeric_sequence(size_t const len = 500,
                               number_type const min = std::numeric_limits<number_type>::lowest(),
                               number_type const max = std::numeric_limits<number_type>::max(),
                               size_t const seed = 0)
{
    std::mt19937_64 engine(seed);
    std::uniform_int_distribution<size_t> dist{min, max};

    auto gen = [&dist, &engine]()
    {
        return dist(engine);
    };
    std::vector<number_type> sequence(len);
    std::ranges::generate(sequence, gen);

    return sequence;
}

template <typename alphabet_t>
auto generate_sequence_pairs(size_t const size, size_t const set_size, size_t const size_variance = 0)
{
    using sequence_t = std::vector<alphabet_t>;

    seqan3::test::random_sequence_generator<sequence_t> generator{size, size_variance};
    std::mt19937_64 random_engine{0};

    std::vector<std::pair<sequence_t, sequence_t>> sequence_pairs(set_size);
    std::ranges::generate(sequence_pairs,
                          [&]() -> std::pair<sequence_t, sequence_t>
                          {
                              return {generator(random_engine), generator(random_engine)};
                          });
    return sequence_pairs;
}

#ifdef SEQAN3_HAS_SEQAN2
template <typename alphabet_t>
auto generate_sequence_seqan2(size_t const len = 500, size_t const variance = 0, size_t const seed = 0)
{
    std::mt19937 gen(seed);
    std::uniform_int_distribution<size_t> dis_alpha(0, seqan2::ValueSize<alphabet_t>::VALUE - 1);
    std::uniform_int_distribution<size_t> dis_length(len - variance, len + variance);

    seqan2::String<alphabet_t> sequence;
    size_t length = dis_length(gen);

    for (size_t l = 0; l < length; ++l)
        appendValue(sequence, static_cast<alphabet_t>(dis_alpha(gen)));

    return sequence;
}

template <typename alphabet_t>
auto generate_sequence_pairs_seqan2(size_t const sequence_length,
                                    size_t const set_size,
                                    size_t const sequence_variance = 0)
{
    using sequence_t = decltype(generate_sequence_seqan2<alphabet_t>());

    seqan2::StringSet<sequence_t> vec1;
    seqan2::StringSet<sequence_t> vec2;

    for (unsigned i = 0; i < set_size; ++i)
    {
        sequence_t seq1 = generate_sequence_seqan2<alphabet_t>(sequence_length, sequence_variance, i);
        sequence_t seq2 = generate_sequence_seqan2<alphabet_t>(sequence_length, sequence_variance, i + set_size);
        appendValue(vec1, seq1);
        appendValue(vec2, seq2);
    }

    return std::tuple{vec1, vec2};
}
#endif // generate seqan2 data.

} // namespace seqan3::test
