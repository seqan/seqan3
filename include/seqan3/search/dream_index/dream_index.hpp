// -----------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2019, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2019, Knut Reinert & MPI für molekulare Genetik
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at: https://github.com/seqan/seqan3/blob/master/LICENSE
// -----------------------------------------------------------------------------------------------------

/*!\file
 * \author Enrico Seiler <enrico.seiler AT fu-berlin.de>
 * \brief Provides seqan3::dream_index.
 */

#pragma once

#include <seqan3/alphabet/nucleotide/dna4.hpp>
#include <seqan3/range/view/kmer_hash.hpp>
#include <seqan3/range/view/persist.hpp>
#include <seqan3/search/dream_index/binning_directory.hpp>
#include <seqan3/search/dream_index/concept.hpp>
#include <seqan3/std/ranges>

namespace seqan3
{

//!\brief The default DREAM index configuration.
struct dream_index_default_traits
{
    //!\brief The alphabet of the indexed text.
    using alphabet_t = dna4;
    //!\brief The used bitvector layout.
    using bitvector_strategy = uncompressed;
    //using shape_strategy = normal; // Tags  SHould be shape later
    //!\brief The used binning data structure.
    using directory_strategy = ibf; // Tags
};

// Bitvector strategy(Un-/Compressed/Chunked), Shape strategy(Normal,Offset,Minimizer,Gapped), Hash strategy(IBF/Direct) => traits

/*!\brief The DREAM index.
 * \tparam dream_index_traits The configuration to use. Must model seqan3::DreamIndexTraits. Defaults to seqan3::dream_index_default_traits.
 * \details \todo Write me.
 */
template <DreamIndexTraits dream_index_traits = dream_index_default_traits>
class dream_index
{
private:
    // using hash_strategy = normal; // Implement via Shapes
    //!\brief The used seqan3::binning_directory specialisation.
    using directory_type = binning_directory<typename dream_index_traits::directory_strategy,
                                             typename dream_index_traits::bitvector_strategy>;
    //!\brief The alphabet type to build the index for.
    using alphabet_t = typename dream_index_traits::alphabet_t;
    //!\brief The number of bins.
    size_t bins;
    //!\brief The k-mer size.
    size_t k; // TODO What about resizing ... window size, offset, etc
    //!\brief The size of the bitvector.
    size_t bits;
    //!\brief The binning_directory.
    directory_type directory;
public:
    /*!\name Constructors, destructor and assignment
     * \{
     */
    dream_index() = default;                                //!< Default constructor.
    dream_index(dream_index const &) = default;             //!< Copy constructor.
    dream_index & operator=(dream_index const &) = default; //!< Move constructor.
    dream_index(dream_index &&) = default;                  //!< Copy assignment.
    dream_index & operator=(dream_index &&) = default;      //!< Move assignment.
    ~dream_index() = default;                               //!< Destructor.

    // TODO Depending on strategy, we need the bits (ibf) or they are fixed (direct)
    /*!\brief Construct a DREAM index given the number of bin, k and the bitvector size.
     * \param b_    The number of bins.
     * \param k_    The k-mer size.
     * \param bits_ The bitvector size.
     */
    dream_index(size_t b_, size_t k_, size_t bits_): bins(b_), k(k_), bits(bits_)
    {
        directory = directory_type(bins, bits);
        // hash_v = hash_type(bins, bits);
    }
    //!\}

    // TODO  this may be a collection?
    /*!\brief Insert a text into a specific bin.
     * \tparam text_t The type of the text. Must model std::ForwardRange and the alphabet type of the text must be equal to alphabet_t.
     * \param bin  The bin to insert the data into.
     * \param text The text to process.
     */
    template <ForwardRange text_t>
    //!\cond
        requires std::Same<innermost_value_type_t<text_t>, alphabet_t>
    //!\endcond
    void insert_data(size_t bin, text_t const & text)
    {
        std::vector<size_t> tmp = text | view::kmer_hash(k);
        for (size_t val : tmp)
        {
            directory.set(val, bin);
        }
    }

    /*!\brief Count the k-mers of a query in all bins.
     * \tparam query_t The type of the query. Must model std::ForwardRange and the alphabet type of the query must be equal to alphabet_t.
     * \param query The query to count the k-mers for.
     * \returns A std::vector<size_t> of size bins where each element is the k-mer count for the respecitve bin.
     */
    template <ForwardRange query_t>
    //!\cond
        requires std::Same<innermost_value_type_t<query_t>, alphabet_t>
    //!\endcond
    std::vector<size_t> count(query_t const & query) const noexcept
    {
        std::vector<size_t> result(bins, 0);
        for (sdsl::bit_vector binning_vector : query | view::kmer_hash(k) | std::view::transform([this](size_t const h)
                                                                            {
                                                                                return directory.get(h);
                                                                            }))
        {
            // TODO SIMD
            size_t bin{0};
            for (size_t batch = 0; batch < ((bins + 63) >> 6)/*directory.bin_width*/; ++batch)
            {
                size_t tmp = binning_vector.get_int(batch * 64);
                if (tmp ^ (1ULL<<63))
                {
                    while (tmp > 0)
                    {
                        uint8_t step = sdsl::bits::lo(tmp);
                        bin += step++;
                        tmp >>= step;
                        ++result[bin++];
                    }
                }
                else
                {
                    ++result[bin + 63];
                }
            }
        }
        return result;
    }

    /*!\brief Determine all bins that contain a query.
     * \tparam query_t The type of the query. Must model std::ForwardRange and the alphabet type of the query must be equal to alphabet_t.
     * \param query The query to determine bin membership for.
     * \param errors The maximum number of allowed errors.
     * \returns A std::vector<size_t> where each element is a bin number.
     */
    template <ForwardRange query_t>
    //!\cond
        requires std::Same<innermost_value_type_t<query_t>, alphabet_t>
    //!\endcond
    std::vector<size_t> get_bins(query_t const & query, uint8_t errors = 0) const noexcept
    {
        std::vector<size_t> result{};
        size_t b{0};
        size_t threshold = (errors + 1) * k > std::ranges::size(query) ? 0 : std::ranges::size(query) - (errors + 1) * k + 1; // 0 or 1?

        auto threshold_filter = std::view::filter([&b, threshold] (auto const c)
        {
            ++b;
            return c >= threshold;
        });

        auto return_bin = std::view::transform([&b] (auto const &)
                          {
                              return b;
                          });

        for (auto e : count(query) | view::persist | threshold_filter | return_bin)
            result.push_back(e);

        return result;
    }
};

} // namespace seqan3
