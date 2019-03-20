// -----------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2019, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2019, Knut Reinert & MPI für molekulare Genetik
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at: https://github.com/seqan/seqan3/blob/master/LICENSE
// -----------------------------------------------------------------------------------------------------

/*!\file
 * \author Enrico Seiler <enrico.seiler AT fu-berlin.de>
 * \brief Provides seqan3::binning_directory.
 */

#pragma once

#include <seqan3/search/dream_index/detail/bitvector.hpp>

// This file takes care of how to process the values of a shape over a text
// Generates the positions in the bitvector and sets them.

// How to give any bitvector_strategy?

namespace seqan3
{

//!\brief Tag for direct addressing.
struct direct {};
//!\brief Tag for Interleaved Bloom Filter.
struct ibf {};

//!\cond
template <typename strategy, typename bitvector_t>
class binning_directory {};
//!\endcond

/*!\brief The IBF binning directory.
 * \tparam bitvector_t The underlying bitvector layout.
 */
template <typename bitvector_t>
class binning_directory<ibf, bitvector_t>
{
private:
    //!\brief The number of bins.
    size_t bins;
    //!\brief The size of the bitvector.
    size_t bits;
    //!\brief The number of 64 bit integers used to represent bins.
    size_t bin_width;
    //!\brief How big is bins in multiple of 64.
    size_t block_size;
    //!\brief How many blocks fit in the bitvector.
    size_t block_count;
    //!\brief The number of hash functions.
    size_t number_hashes{3};
    //!\brief The bitvector.
    bitvector<bitvector_t> data;
    //!\brief Precalculated values for hashing.
    std::vector<size_t> precalc_values;
    //!\brief Shift value for hashing.
    static constexpr size_t shift{27};
    //!\brief Seed for hashing.
    static constexpr size_t seed{0x90b45d39fb6da1fa};

public:
    /*!\name Constructors, destructor and assignment
     * \{
     */
    binning_directory() = default;                                      //!< Default constructor.
    binning_directory(binning_directory const &) = default;             //!< Copy constructor.
    binning_directory & operator=(binning_directory const &) = default; //!< Move constructor.
    binning_directory(binning_directory &&) = default;                  //!< Copy assignment.
    binning_directory & operator=(binning_directory &&) = default;      //!< Move assignment.
    ~binning_directory() = default;                                               //!< Destructor.

    // Gets a hash value and sets the bits
    // Gets a hash value and returns the binning vector
    // also the data

    /*!\brief Construct using number of bins and bitvector size.
     * \param bins_ The number of bins.
     * \param bits_ The bitvector size.
     */
    binning_directory(size_t bins_, size_t bits_): bins(bins_), bits(bits_)
    {
        data = bitvector<bitvector_t>(bits);
        bin_width = (bins + 63) >> 6;
        block_size = bin_width << 6;
        block_count = bits / block_size;

        precalc_values.resize(number_hashes);
        for(size_t i = 0; i < number_hashes ; ++i)
            precalc_values[i] = i ^ (15 * seed);
    }
    //!\}

    /*!\brief Perturbs a value and fits it into the vector.
     * \param h The value to process.
     */
    inline constexpr void hash_and_fit(size_t & h) const
    {
        h ^= h >> shift; // Basically Fibonacci hashing
        h %= block_count;
        h *= block_size;
    }

    /*!\brief Inserts a value into a specific bin.
     * \param h   The raw hash value to process.
     * \param bin The bin to insert into.
     */
    inline void set(size_t h, size_t bin)
    {
        // Sets respective positions in bitvector TODO SIMD?
        std::ranges::for_each(precalc_values, [](size_t const & val) {
            size_t idx = val * h;
            hash_and_fit(idx);
            idx += bin;
            data[idx] = 1;
        });
    }

    // Do we want this?
    //!\cond
    inline void set(size_t h, std::vector<size_t> bin)
    {
    }
    //!\endcond

    /*!\brief Determines set membership of a given value.
     * \param h   The raw hash value to process.
     * \returns A sdsl::bit_vector of size bins where each position indicates the bin membership of the value.
     */
    inline sdsl::bit_vector get(size_t h) const
    {
        // returns binning vector
        std::vector<size_t> indices = precalc_values; // TODO this is copied every time (assignment_time = copy_time?)
        sdsl::bit_vector result(bins);

        std::ranges::for_each(indices, [](size_t & idx) { idx *= h; hash_and_fit(h); })

        // SIMD?!! std::span?
        for (size_t batch = 0; batch < bin_width; ++batch)
        {
           size_t tmp = data.get_int(indices[0]);
           indices[0] += 64;

           for(size_t i = 1; i < number_hashes; ++i)
           {
               tmp &= data.get_int(indices[i]);
               indices[i] += 64;
           }

           result.set_int(64 * batch, tmp);
        }
        return result;
    }
};

} // namespace seqan3
