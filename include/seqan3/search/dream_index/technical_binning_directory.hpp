// -----------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2020, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2020, Knut Reinert & MPI für molekulare Genetik
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at: https://github.com/seqan/seqan3/blob/master/LICENSE.md
// -----------------------------------------------------------------------------------------------------

/*!\file
 * \author Enrico Seiler <enrico.seiler AT fu-berlin.de>
 * \brief Provides seqan3::technical_binning_directory.
 */

#pragma once

#include <seqan3/alphabet/nucleotide/dna4.hpp>
#include <seqan3/core/algorithm/detail/execution_handler_parallel.hpp>
#include <seqan3/utility/tuple/concept.hpp>
#include <seqan3/range/views/chunk.hpp>
#include <seqan3/range/views/kmer_hash.hpp>
#include <seqan3/range/views/zip.hpp>
#include <seqan3/search/dream_index/interleaved_bloom_filter.hpp>

#include <seqan3/range/views/async_input_buffer.hpp>

namespace seqan3
{

/*!\addtogroup submodule_dream_index
 * \{
 */
//!\brief Stores paramters to construct the seqan3::interleaved_bloom_filter with.
struct ibf_config
{
    bin_count number_of_bins; //!< The number of bins.
    bin_size size_of_bin; //!< The size of each individual bin.
    hash_function_count number_of_hash_functions; //!< The number of hash functions.
    uint8_t threads{1u}; //!< The number of threads to use functions.
};

enum class hash_variant : uint8_t
{
    kmer,
    minimiser
};

template <auto t, typename hasher_t>
struct hash_proxy;

template <typename hasher_t>
struct hash_proxy<hash_variant::kmer, hasher_t>
{
    size_t kmer_size;
    hasher_t hasher(/*seqan3::views::kmer_hash(seqan3::ungapped{5u})*/); // view is not default constructible

    hash_proxy(size_t kmer_size_) : kmer_size(kmer_size_)/*,
                                    hasher(seqan3::views::kmer_hash(seqan3::ungapped{kmer_size}))*/ {}
};

template <typename hasher_t>
struct hash_proxy<hash_variant::minimiser, hasher_t>
{
    size_t kmer_size;
    size_t window_size;
    hasher_t hasher(/*seqan3::views::minimiser_hash(10u, 5u)*/); // view is not default constructible

    hash_proxy(size_t kmer_size_, size_t window_size_) : kmer_size(kmer_size_), window_size(window_size_)/*,
                                    hasher(seqan3::views::minimiser_hash(kmer_size, window_size))*/ {}
};

/*!\brief The Technical Binning Directory. A data structure that enhances the seqan3::interleaved_bloom_filter by
 * handling sequences as input and query.
 * \tparam data_layout_mode_ Indicates whether the underlying data type is compressed. See seqan3::data_layout.
 * \tparam hash_adaptor_t The type of the view adaptor to generate hash values from a sequence.
 * \tparam alph_t The alphabet type of the technical bins. Must model seqan3::semialphabet.
 * \implements seqan3::cerealisable
 * \sa seqan3::interleaved_bloom_filter
 *
 * \details
 *
 * ### Difference to the Interleaved Bloom Filter
 *
 * In addition to the seqan::interleaved_bloom_filter, the Technical Binning Directory supports construction via
 * a range of sequences.
 * Furthermore, counting k-mer occurrences of a query is supported via the seqan3::counting_agent_type.
 *
 * ### Technical Bins
 *
 * A technical bin is a sequence collection that maps 1-on-1 to the bin in the Technical Binning Directory.
 *
 * ### Querying
 *
 * To query the Technical Binning Directory for a value use the seqan3::membership_agent.
 * To query and count the Technical Binning Directory for all hash values within a sequence, use the
 * seqan3::counting_agent_type.
 *
 * ### Compression
 *
 * The Technical Binning Directory can be compressed by passing `data_layout::compressed` as template argument.
 * The compressed `seqan3::technical_binning_directory<seqan3::data_layout::compressed>` can only be constructed from a
 * `seqan3::technical_binning_directory`, in which case the underlying bitvector is compressed.
 * The compressed Technical Binning Directory is immutable, i.e. only querying is supported.
 *
 * ### Thread safety
 *
 * The Technical Binning Directory promises the basic thread-safety by the STL that all
 * calls to `const` member functions are safe from multiple threads (as long as no thread calls
 * a non-`const` member function at the same time).
 *
 * Additionally, concurrent calls to `emplace` are safe iff each thread handles a multiple of wordsize (=64) many bins.
 * For example, calls to `emplace` from multiple threads are safe if `thread_1` accesses bins 0-63, `thread_2` bins
 * 64-127, and so on.
 */
template <data_layout data_layout_mode_ = data_layout::uncompressed,
          typename hash_adaptor_t_ = decltype(seqan3::views::kmer_hash(seqan3::ungapped(5u))),
          semialphabet alph_t = dna4>
class technical_binning_directory : public interleaved_bloom_filter<data_layout_mode_>
{
private:
    //!\cond
    template <data_layout data_layout_mode, typename friend_view_t, semialphabet friend_alph_t>
        requires std::same_as<friend_view_t, hash_adaptor_t_>
    friend class technical_binning_directory;

    template <std::integral value_t>
    class counting_agent_type;
    //!\endcond

    //!\brief The type of the underlying IBF.
    using base_t = interleaved_bloom_filter<data_layout_mode_>;
    //!\brief The adaptor to use for generating hash values from a sequence.
    hash_adaptor_t_ hash_adaptor;

public:
    //!\brief Indicates whether the Technical Binning Directory is compressed.
    static constexpr data_layout data_layout_mode = data_layout_mode_;
    //!\brief The type of the hash adaptor.
    using hash_adaptor_t = hash_adaptor_t_;

    technical_binning_directory(): hash_adaptor(seqan3::views::kmer_hash(seqan3::ungapped{5u})) {};
    technical_binning_directory(technical_binning_directory const &) = default; //!< Defaulted.
    technical_binning_directory & operator=(technical_binning_directory const &) = default; //!< Defaulted.
    technical_binning_directory(technical_binning_directory &&) = default; //!< Defaulted.
    technical_binning_directory & operator=(technical_binning_directory &&) = default; //!< Defaulted.
    ~technical_binning_directory() = default; //!< Defaulted.


    /*!\brief Construct an uncompressed Technical Binning Directory.
     * \tparam rng_t The type of the technical bins.
     * \param technical_bins The technical bins.
     * \param hash_adaptor A hash adaptor determining how to hash the input sequences.
     * \param cfg A seqan3::ibf_config.
     *
     * \attention This constructor can only be used to construct **uncompressed** Technical Binning Directories.
     *
     * \details
     *
     * ### Example
     *
     * \include test/snippet/search/dream_index/technical_binning_directory.cpp
     */
    template <std::ranges::range rng_t>
    //!\cond
        requires (data_layout_mode == data_layout::uncompressed)
    //!\endcond
    technical_binning_directory(rng_t && technical_bins,
                                hash_adaptor_t hash_adaptor,
                                ibf_config const & cfg)
        : base_t(cfg.number_of_bins, cfg.size_of_bin, cfg.number_of_hash_functions),
          hash_adaptor(std::move(hash_adaptor))
    {
        static_assert(range_dimension_v<rng_t> == 2 || range_dimension_v<rng_t> == 3,
                      "Technical bins must be given as range of ranges or range of ranges of ranges (!).");
        static_assert(std::ranges::input_range<rng_t>, "Technical bins must model input_range.");
        static_assert(std::ranges::input_range<std::ranges::range_reference_t<rng_t>>,
                      "Individual bins must model input_range.");
        // static_assert(semialphabet<range_innermost_value_t<rng_t>>, "The content of a bin must model semialphabet.");
        // static_assert(semialphabet<range_innermost_value_t<rng_t::traits_type::sequence_type>>);

        size_t const number_of_bins = cfg.number_of_bins.get();

        using rng_difference_t = std::ranges::range_difference_t<rng_t>;
        if (std::ranges::distance(technical_bins) > static_cast<rng_difference_t>(number_of_bins))
            throw std::logic_error("There are more bins in the input data set than "
                                   "the binning directory is configured to handle.");

        auto worker = [&] (auto && zipped_view, auto &&)
        {
            auto hash_adaptor_copy = hash_adaptor;
            for (auto && [technical_bin, bin_number] : zipped_view)
            {
                bin_index const idx{bin_number};
                if constexpr(seqan3::tuple_like<std::ranges::range_reference_t<std::ranges::range_reference_t<rng_t>>>)
                {
                    for (auto && [seq] : technical_bin)
                    {
                        for (auto && hash : seq | hash_adaptor_copy)
                            this->emplace(hash, idx);
                    }
                }
                else
                {
                    for (auto && hash : technical_bin | hash_adaptor_copy)
                        this->emplace(hash, idx);
                }
            }
        };

        auto worker_async = [&] (auto && zipped_view, auto &&)
        {
            auto hash_adaptor_copy = hash_adaptor;
            for (auto && [technical_bin, bin_number] : zipped_view)
            {
                bin_index const idx{bin_number};
                if constexpr(seqan3::tuple_like<std::ranges::range_reference_t<std::ranges::range_reference_t<rng_t>>>)
                {
                    for (auto && [seq] : technical_bin | seqan3::views::async_input_buffer(2))
                    {
                        for (auto && hash : seq | hash_adaptor_copy)
                            this->emplace(hash, idx);
                    }
                }
                else
                {
                    for (auto && hash : technical_bin | hash_adaptor_copy)
                        this->emplace(hash, idx);
                }
            }
        };

        // A single thread may handle between 8 and 64 bins.
        size_t const chunk_size = std::clamp<size_t>(std::bit_ceil(number_of_bins / cfg.threads),
                                                     8u,
                                                     64u);
        // If at least half of the threads are idle, we will make the I/O asynchronous for them.
        size_t const threads = chunk_size * cfg.threads <= number_of_bins ? cfg.threads :
                                   chunk_size * cfg.threads / 2 <= number_of_bins ? cfg.threads : cfg.threads / 2;
        auto chunked_view = views::zip(technical_bins, std::views::iota(0u)) | views::chunk(chunk_size);
        detail::execution_handler_parallel executioner{threads};
        if (threads == cfg.threads / 2)
            executioner.bulk_execute(worker_async, std::move(chunked_view), [](){});
        else
            executioner.bulk_execute(worker, std::move(chunked_view), [](){});
    }

    //!\cond
    // Constructor for IBF App...
    // I need to construct an empty TBD for serialisation.
    // I need to provide the config (do I?) and the hash_adaptor(not serialisable), but can ignore the technical_bins.
    template <std::ranges::range rng_t>
        requires (data_layout_mode == data_layout::uncompressed)
    technical_binning_directory(rng_t && technical_bins,
                                hash_adaptor_t hash_adaptor,
                                ibf_config const & cfg,
                                bool)
        : base_t(cfg.number_of_bins, cfg.size_of_bin, cfg.number_of_hash_functions),
          hash_adaptor(std::move(hash_adaptor))
    {
        (void) technical_bins;
    }
    //!\endcond

    /*!\brief Construct a compressed Technical Binning Directory.
     * \param[in] tbd The uncompressed seqan3::technical_binning_directory.
     *
     * \attention This constructor can only be used to construct **compressed** Technical Binning Directorys.
     *
     * \details
     *
     * ### Example
     *
     * \include test/snippet/search/dream_index/technical_binning_directory.cpp
     */
    technical_binning_directory(technical_binning_directory<data_layout::uncompressed,
                                                            hash_adaptor_t,
                                                            alph_t> && tbd)
    //!\cond
        requires (data_layout_mode == data_layout::compressed)
    //!\endcond
        : base_t(tbd), hash_adaptor(tbd.hash_adaptor)
    {}

    /*!\name Lookup
     * \{
     */
    /*!\brief Returns seqan3::technical_binning_directory::counting_agent_type to be used for lookup.
     * \attention Calling seqan3::technical_binning_directory::increase_bin_number_to invalidates all
     * seqan3::technical_binning_directory::counting_agent_type constructed for this Technical Binning Directory.
     *
     * \details
     *
     * ### Example
     *
     * \include test/snippet/search/dream_index/counting_agent.cpp
     * \sa seqan3::technical_binning_directory::counting_agent_type::bulk_contains
     */
    template <std::integral value_t = size_t>
    counting_agent_type<value_t> counting_agent() const
    {
        return counting_agent_type<value_t>{*this};
    }
    //!\}
};

/*!\name Type deduction guides
 * \{
 */
//!\brief Deduces the template.
template <typename rng_t, typename hash_adaptor_t>
technical_binning_directory(rng_t, hash_adaptor_t, ibf_config) ->
        technical_binning_directory<data_layout::uncompressed, hash_adaptor_t, range_innermost_value_t<rng_t>>;
//!\}

/*!\brief Manages counting queries for the seqan3::technical_binning_directory.
 * \attention Calling seqan3::technical_binning_directory::increase_bin_number_to invalidates the counting_agent_type.
 *
 * \details
 *
 * ### Example
 *
 * \include test/snippet/search/dream_index/counting_agent.cpp
 */
template <data_layout data_layout_mode, typename hash_adaptor_t, semialphabet alph_t>
template <std::integral value_t>
class technical_binning_directory<data_layout_mode, hash_adaptor_t, alph_t>::counting_agent_type
{
private:
    //!\brief The type of the augmented seqan3::technical_binning_directory.
    using tbd_t = technical_binning_directory<data_layout_mode, hash_adaptor_t, alph_t>;

    //!\brief A pointer to the augmented seqan3::technical_binning_directory.
    tbd_t const * tbd_ptr;

    //!\brief Store a seqan3::technical_binning_directory::membership_agent to call `bulk_contains`.
    decltype(std::declval<tbd_t>().membership_agent()) membership_agent;

public:
    /*!\name Constructors, destructor and assignment
     * \{
     */
    counting_agent_type() = default; //!< Defaulted.
    counting_agent_type(counting_agent_type const &) = default; //!< Defaulted.
    counting_agent_type & operator=(counting_agent_type const &) = default; //!< Defaulted.
    counting_agent_type(counting_agent_type &&) = default; //!< Defaulted.
    counting_agent_type & operator=(counting_agent_type &&) = default; //!< Defaulted.
    ~counting_agent_type() = default; //!< Defaulted.

    /*!\brief Construct a counting_agent_type for an existing seqan3::technical_binning_directory.
     * \private
     * \param tbd The seqan3::technical_binning_directory.
     */
    counting_agent_type(tbd_t const & tbd) : tbd_ptr(std::addressof(tbd)), membership_agent(tbd)
    {
        result_buffer.resize(tbd_ptr->bin_count());
    };
    //!\}

    //!\brief Stores the result of bulk_contains().
    counting_vector<value_t> result_buffer;

    /*!\name Counting
     * \{
     */
    /*!\brief Determines set membership of a given query.
     * \param[in] query The sequence to process.
     *
     * \attention The result of this function must always be bound via reference, e.g. `auto &` to prevent copying.
     * \attention Sequential calls to this function invalidate the previously returned reference.
     *
     * \details
     *
     * ### Example
     *
     * \include test/snippet/search/dream_index/counting_agent.cpp
     *
     * ### Thread safety
     *
     * Concurrent invocations of this function are not thread safe, please create a seqan3::counting_agent_type for each thread.
     */
    template <std::ranges::range query_t>
    counting_vector<value_t> const & count_query(query_t const & query) & noexcept
    {
        assert(tbd_ptr != nullptr);
        assert(result_buffer.size() == tbd_ptr->bin_count());

        static_assert(std::ranges::input_range<query_t>, "The query must model input_range.");
        static_assert(std::convertible_to<range_innermost_value_t<query_t>, alph_t>,
                      "The alphabet of the query must be convertible to the alphabet of the Technical Binning "
                      "Directory.");

        std::ranges::fill(result_buffer, 0);

        for (auto && hash : query | tbd_ptr->hash_adaptor)
            result_buffer += membership_agent.bulk_contains(hash);

        return result_buffer;
    }

    /*!\brief Determines set membership of a given query and counts number of searched hashes.
     * \param[in] query The sequence to process.
     *
     * \attention The result of this function must always be bound via reference, e.g. `auto &` to prevent copying.
     * \attention Sequential calls to this function invalidate the previously returned reference.
     *
     * \details
     *
     * ### Example
     *
     * \include test/snippet/search/dream_index/counting_agent.cpp
     *
     * ### Thread safety
     *
     * Concurrent invocations of this function are not thread safe, please create a seqan3::counting_agent_type for each thread.
     */
    template <std::ranges::range query_t>
    std::pair<counting_vector<value_t> const &, value_t const> count_query(query_t const & query, bool) & noexcept // TODO no second parameter please
    {
        assert(tbd_ptr != nullptr);
        assert(result_buffer.size() == tbd_ptr->bin_count());

        static_assert(std::ranges::input_range<query_t>, "The query must model input_range.");
        static_assert(std::convertible_to<range_innermost_value_t<query_t>, alph_t>,
                      "The alphabet of the query must be convertible to the alphabet of the Technical Binning "
                      "Directory.");

        std::ranges::fill(result_buffer, 0);

        value_t count{};
        for (auto && hash : query | tbd_ptr->hash_adaptor)
        {
            result_buffer += membership_agent.bulk_contains(hash);
            ++count;
        }

        return {result_buffer, count};
    }

    // `bulk_contains` cannot be called on a temporary, since the object the returned reference points to
    // is immediately destroyed.
    counting_vector<value_t> const & count_query(size_t const value) && noexcept = delete;
    //!\}

};

//!\}

} // namespace seqan3
