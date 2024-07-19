// SPDX-FileCopyrightText: 2006-2024 Knut Reinert & Freie Universität Berlin
// SPDX-FileCopyrightText: 2016-2024 Knut Reinert & MPI für molekulare Genetik
// SPDX-License-Identifier: BSD-3-Clause

/*!\file
 * \author Tobias Loka <tobias.loka AT hpi.de>
 * \brief Provides seqan3::bloom_filter.
 */

#pragma once

#include <seqan3/search/dream_index/interleaved_bloom_filter.hpp>

namespace seqan3
{

/*!\brief The Bloom Filter. A data structure that efficiently answers set-membership queries.
 * \tparam data_layout_mode_ Indicates whether the underlying data type is compressed. See seqan3::data_layout.
 * \implements seqan3::cerealisable
 * \ingroup utility_bloom_filter
 *
 * \details
 *
 * ### Bloom Filter (BF)
 *
 * The [Bloom Filter](https://en.wikipedia.org/wiki/Bloom_filter) is a probabilistic data structure.
 * A Bloom Filter can be thought of as a bitvector of length `n` and `h` hash functions and is used to determine set
 * membership. To insert data, the data is hashed by the `h` hash functions (returning values in `[0, n)`) and the
 * corresponding `h` positions in the bitvector are set to `1`. To query data, i.e. to determine whether the query
 * belongs to the set the Bloom Filter was built for, the query is hashed by the same `h` hash functions and the
 * corresponding positions are checked. If all `h` positions contain a `1`, the query is (probably) in the data set.
 * Since the Bloom Filter has variable length, the hashing is not bijective, i.e. it may return true for a set
 * membership query even though the query was never inserted into the Bloom Filter. Note that the Bloom Filter
 * will always return `true` if the query was inserted, i.e. there may be false positives, but no false negatives.
 *
 * ### Querying
 *
 * To query the Bloom Filter for a value, call `seqan3::bloom_filter::contains` which returns
 * true if the k-mer hash is present in the index, and false if the hash is not present.
 * The value is a hash value of the k-mer to check membership for.
 *
 * To query the Bloom Filter for a range of values, call `seqan3::bloom_filter::count` which returns the
 * number of k-mer hits in the Bloom Filter for the given range of values.
 *
 * Please note the results are based on a heuristic data structure and with a certain probability (depending on
 * the selected size of the bit vector) you may receive a false positive result.
 *
 * ### Differences to the Interleaved Bloom Filter (IBF)
 *
 * While the Bloom Filter provides a single linear bit vector to represent the underlying data, the Interleaved Bloom
 * Filter provides a data structure that combines a set of Bloom Filters to enable efficient queries to multiple
 * fractions of the data. In doing so, the Interleaved Bloom Filter can not only answer whether a hash
 * value is present in the data, but also provides information in which fraction of the data it occurs.
 * The design of the Interleaved Bloom Filter is particularly useful when the underlying data is systematically
 * structured; for example, if each fraction of the data represents a specific set of organisms. Important
 * applications of the Interleaved Bloom Filter include taxonomic classification of sequencing data, or prefiltering
 * of specific fractions of an input data set to enable more efficient in-depth analysis.
 * The Bloom Filter, on the other hand, is useful if the database does not contain any underlying structure, or it is
 * not relevant for the analysis. A typical application is the removal of host sequences or different types of
 * contamination where it is usually not of interest which part of the database was matched. In such cases, the Bloom
 * Filter provides a lighter data structure and a more simple interface (for example, the use of agents for determining
 * and counting membership is not necessary in this case).
 *
 * ### Compression
 *
 * The Bloom Filter can be compressed by passing `seqan3::data_layout::compressed` as template argument.
 * The compressed `seqan3::bloom_filter<seqan3::data_layout::compressed>` can only be constructed from a
 * `seqan3::bloom_filter`, in which case the underlying bitvector is compressed.
 * The compressed Bloom Filter is immutable, i.e. only querying is supported.
 *
 * ### Thread safety
 *
 * The Bloom Filter promises the basic thread-safety by the STL that all
 * calls to `const` member functions are safe from multiple threads (as long as no thread calls
 * a non-`const` member function at the same time).
 *
 * \sa seqan3::interleaved_bloom_filter
 *
 */
template <data_layout data_layout_mode_ = data_layout::uncompressed>
class bloom_filter
{
private:
    //!\cond
    template <data_layout data_layout_mode>
    friend class bloom_filter;
    //!\endcond

    //!\brief The underlying datatype to use.
    using data_type =
        std::conditional_t<data_layout_mode_ == data_layout::uncompressed, sdsl::bit_vector, sdsl::sd_vector<>>;

    //!\brief The size of the underlying bit vector in bits.
    size_t size_in_bits{};
    //!\brief The number of bits to shift the hash value before doing multiplicative hashing.
    size_t hash_shift{};
    //!\brief The number of hash functions.
    size_t hash_funs{};
    //!\brief The bitvector.
    data_type data{};
    //!\brief Precalculated seeds for multiplicative hashing. We use large irrational numbers for a uniform hashing.
    static constexpr std::array<size_t, 5> hash_seeds{13'572'355'802'537'770'549ULL, // 2**64 / (e/2)
                                                      13'043'817'825'332'782'213ULL, // 2**64 / sqrt(2)
                                                      10'650'232'656'628'343'401ULL, // 2**64 / sqrt(3)
                                                      16'499'269'484'942'379'435ULL, // 2**64 / (sqrt(5)/2)
                                                      4'893'150'838'803'335'377ULL}; // 2**64 / (3*pi/5)

    /*!\brief Perturbs a value and fits it into the vector.
     * \param h The value to process.
     * \param seed The seed to use.
     * \returns A hashed value representing a position within the bounds of `data`.
     * \sa https://probablydance.com/2018/06/16/
     * \sa https://lemire.me/blog/2016/06/27
     */
    inline constexpr size_t hash_and_fit(size_t h, size_t const seed) const
    {
        h *= seed;
        h ^= h >> hash_shift;               // XOR and shift higher bits into lower bits
        h *= 11'400'714'819'323'198'485ULL; // = 2^64 / golden_ration, to expand h to 64 bit range
                                            // Use fastrange (integer modulo without division) if possible.
#ifdef __SIZEOF_INT128__
        h = static_cast<uint64_t>((static_cast<__uint128_t>(h) * static_cast<__uint128_t>(size_in_bits)) >> 64);
#else
        h %= size_in_bits;
#endif
        return h;
    }

public:
    //!\brief Indicates whether the Bloom Filter is compressed.
    static constexpr data_layout data_layout_mode = data_layout_mode_;

    /*!\name Constructors, destructor and assignment
     * \{
     */
    bloom_filter() = default;                                 //!< Defaulted.
    bloom_filter(bloom_filter const &) = default;             //!< Defaulted.
    bloom_filter & operator=(bloom_filter const &) = default; //!< Defaulted.
    bloom_filter(bloom_filter &&) = default;                  //!< Defaulted.
    bloom_filter & operator=(bloom_filter &&) = default;      //!< Defaulted.
    ~bloom_filter() = default;                                //!< Defaulted.

    /*!\brief Construct an uncompressed Bloom Filter.
     * \param size The bit vector size (in bits).
     * \param funs The number of hash functions. Default 2. At least 1, at most 5.
     *
     * \attention This constructor can only be used to construct **uncompressed** Bloom Filters.
     *
     * \details
     *
     * ### Example
     *
     * \include test/snippet/utility/bloom_filter/bloom_filter_constructor.cpp
     */
    bloom_filter(seqan3::bin_size size, seqan3::hash_function_count funs = seqan3::hash_function_count{2u})
        requires (data_layout_mode == data_layout::uncompressed)
    {
        size_in_bits = size.get();
        hash_funs = funs.get();

        if (hash_funs == 0 || hash_funs > 5)
            throw std::logic_error{"The number of hash functions must be > 0 and <= 5."};
        if (size_in_bits == 0)
            throw std::logic_error{"The size of a bloom filter must be > 0."};

        hash_shift = std::countl_zero(size_in_bits);
        data = sdsl::bit_vector(size_in_bits);
    }

    /*!\brief Construct a compressed Bloom Filter.
     * \param[in] bf The uncompressed seqan3::bloom_filter.
     *
     * \attention This constructor can only be used to construct **compressed** Bloom Filters.
     *
     * \details
     *
     * ### Example
     *
     * \include test/snippet/utility/bloom_filter/bloom_filter_constructor_compressed.cpp
     */
    bloom_filter(bloom_filter<data_layout::uncompressed> const & bf)
        requires (data_layout_mode == data_layout::compressed)
    {
        std::tie(size_in_bits, hash_shift, hash_funs) = std::tie(bf.size_in_bits, bf.hash_shift, bf.hash_funs);

        data = sdsl::sd_vector<>{bf.data};
    }
    //!\}

    /*!\name Modifiers
     * \{
     */
    /*!\brief Inserts a value into the Bloom Filter.
     * \param[in] value The raw numeric value to process.
     *
     * \attention This function is only available for **uncompressed** Bloom Filters.
     *
     * \details
     *
     * ### Example
     *
     * \include test/snippet/utility/bloom_filter/bloom_filter_emplace.cpp
     */
    void emplace(size_t const value) noexcept
        requires (data_layout_mode == data_layout::uncompressed)
    {
        for (size_t i = 0; i < hash_funs; ++i)
        {
            size_t idx = hash_and_fit(value, hash_seeds[i]);
            assert(idx < data.size());
            data[idx] = 1;
        };
    }

    /*!\brief Remove all values from the Bloom Filter by setting all bits to 0.
     *
     * \attention This function is only available for **uncompressed** Bloom Filters.
     *
     * \details
     *
     * While all values are removed from the vector, the size of the Bloom Filter is not changed.
     *
     * ### Example
     *
     * \include test/snippet/utility/bloom_filter/bloom_filter_reset.cpp
     */
    void reset() noexcept
        requires (data_layout_mode == data_layout::uncompressed)
    {
        sdsl::util::_set_zero_bits(data);
    }
    //!\}

    /*!\name Lookup
     * \{
     */
    /*!\brief Check whether a value is present in the Bloom Filter.
     * \param[in] value The raw numeric value to process.
     *
     * \attention This function is only available for **uncompressed** Bloom Filters.
     *
     * \details
     *
     * ### Example
     *
     * \include test/snippet/utility/bloom_filter/bloom_filter_contains.cpp
     */
    bool contains(size_t const value) const noexcept
    {
        for (size_t i = 0; i < hash_funs; i++)
        {
            size_t idx = hash_and_fit(value, hash_seeds[i]);
            assert(idx < data.size());
            if (data[idx] == 0)
                return false;
        }
        return true;
    }
    //!\}

    /*!\name Counting
     * \{
     */
    /*!\brief Counts the occurrences for all values in a range.
     * \tparam value_range_t The type of the range of values. Must model std::ranges::input_range. The reference type
     *                       must model std::unsigned_integral.
     * \param[in] values The range of values to process.
     *
     * \details
     *
     * ### Example
     *
     * \include test/snippet/utility/bloom_filter/bloom_filter_count.cpp
     *
     * ### Thread safety
     *
     * Concurrent invocations of this function are thread safe.
     */
    template <std::ranges::range value_range_t>
    size_t count(value_range_t && values) const noexcept
    {
        static_assert(std::ranges::input_range<value_range_t>, "The values must model input_range.");
        static_assert(std::unsigned_integral<std::ranges::range_value_t<value_range_t>>,
                      "An individual value must be an unsigned integral.");

        size_t result = 0;

        for (auto && value : values)
            result += contains(value);

        return result;
    }
    //!\}

    /*!\name Capacity
     * \{
     */
    /*!\brief Returns the number of hash functions used in the Bloom Filter.
     * \returns The number of hash functions.
     */
    size_t hash_function_count() const noexcept
    {
        return hash_funs;
    }

    /*!\brief Returns the size of the underlying bitvector.
     * \returns The size in bits of the underlying bitvector.
     */
    size_t bit_size() const noexcept
    {
        return size_in_bits;
    }
    //!\}

    /*!\name Comparison operators
     * \{
     */
    /*!\brief Test for equality.
     * \param[in] lhs A `seqan3::bloom_filter`.
     * \param[in] rhs `seqan3::bloom_filter` to compare to.
     * \returns `true` if equal, `false` otherwise.
     */
    friend bool operator==(bloom_filter const & lhs, bloom_filter const & rhs) noexcept
    {
        return std::tie(lhs.size_in_bits, lhs.hash_shift, lhs.hash_funs, lhs.data)
            == std::tie(rhs.size_in_bits, rhs.hash_shift, rhs.hash_funs, rhs.data);
    }

    /*!\brief Test for inequality.
     * \param[in] lhs A `seqan3::bloom_filter`.
     * \param[in] rhs `seqan3::bloom_filter` to compare to.
     * \returns `true` if unequal, `false` otherwise.
     */
    friend bool operator!=(bloom_filter const & lhs, bloom_filter const & rhs) noexcept
    {
        return !(lhs == rhs);
    }
    //!\}

    /*!\name Access
     * \{
     */
    /*!\brief Provides direct, unsafe access to the underlying data structure.
     * \returns A reference to an SDSL bitvector.
     *
     * \details
     *
     * \noapi{The exact representation of the data is implementation defined.}
     */
    constexpr data_type & raw_data() noexcept
    {
        return data;
    }

    //!\copydoc raw_data()
    constexpr data_type const & raw_data() const noexcept
    {
        return data;
    }
    //!\}

    /*!\cond DEV
     * \brief Serialisation support function.
     * \tparam archive_t Type of `archive`; must satisfy seqan3::cereal_archive.
     * \param[in] archive The archive being serialised from/to.
     *
     * \attention These functions are never called directly, see \ref serialisation for more details.
     */
    template <cereal_archive archive_t>
    void CEREAL_SERIALIZE_FUNCTION_NAME(archive_t & archive)
    {
        archive(size_in_bits);
        archive(hash_shift);
        archive(hash_funs);
        archive(data);
    }
    //!\endcond
};

} // namespace seqan3
