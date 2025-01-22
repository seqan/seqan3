// SPDX-FileCopyrightText: 2006-2025 Knut Reinert & Freie Universität Berlin
// SPDX-FileCopyrightText: 2016-2025 Knut Reinert & MPI für molekulare Genetik
// SPDX-License-Identifier: BSD-3-Clause
/*!\file
 * \author Mitra Darvish <mitra.darvish AT fu-berlin.de>
 * \brief Provides a SeqAn2 minimiser hash.
 */
#pragma once

#include <seqan3/core/platform.hpp>

#ifdef SEQAN3_HAS_SEQAN2
#include <seqan/seq_io.h>

//!\brief Strong type for passing the window size.
struct window { uint64_t v; };
//!\brief Strong type for passing the kmer size.
struct kmer { uint64_t v; };
//!\brief Strong type for passing number of bins.
struct bins { uint64_t v; };
//!\brief Strong type for passing number of bits.
struct bits { uint64_t v; };
//!\brief Strong type for passing number of hash functions.
struct hashes { uint64_t v; };

template<typename shape_t>
class minimiser
{
private:
    //!\brief The alphabet type.
    using alphabet_t = seqan2::Dna;
    //!\brief The text type.
    using text_t = seqan2::String<alphabet_t>;
    //!\brief The type of the complemented text.
    using complement_t = seqan2::ModifiedString<text_t const, seqan2::ModComplementDna>;
    //!\brief The type of the reverse complemented text.
    using reverse_complement_t = seqan2::ModifiedString<complement_t, seqan2::ModReverse>;

    //!\brief The window size of the minimiser.
    uint64_t w{};
    //!\brief The size of the k-mers.
    uint64_t k{};
    //!\brief Random but fixed value to xor k-mers with. Counteracts consecutive minimisers.
    uint64_t seed{};

    //!\brief Shape for computing the forward strand k-mers.
    shape_t forward_shape{};
    //!\brief Shape for computing the reverse strand k-mers.
    shape_t reverse_shape{};

    //!\brief Stores the k-mer hashes of the forward strand.
    std::vector<uint64_t> forward_hashes;
    //!\brief Stores the k-mer hashes of the reverse complement strand.
    std::vector<uint64_t> reverse_hashes;

public:

    //!\brief Stores the hashes of the minimisers.
    std::vector<uint64_t> minimiser_hash;
    //!\brief Stores the begin positions of the minimisers.
    std::vector<uint64_t> minimiser_begin;
    //!\brief Stores the end positions of the minimisers.
    std::vector<uint64_t> minimiser_end;

    minimiser() = default;                                        //!< Defaulted
    minimiser(minimiser const &) = default;                       //!< Defaulted
    minimiser(minimiser &&) = default;                            //!< Defaulted
    minimiser & operator=(minimiser const &) = default;           //!< Defaulted
    minimiser & operator=(minimiser &&) = default;                //!< Defaulted
    ~minimiser() = default;                                       //!< Defaulted

    /*!\brief Constructs a minimiser from given k-mer, window size and a seed.
     * \param[in] w_     The window size.
     * \param[in] k_     The k-mer size.
     * \param[in] shape_ The shape to use.
     * \param[in] seed_  The seed to use. Default: 0x8F3F73B5CF1C9ADE.
     */
    minimiser(window const w_, kmer const k_, shape_t shape_, uint64_t const seed_ = 0x8F3F73B5CF1C9ADE) :
        w{w_.v}, k{k_.v}, seed{seed_}, forward_shape{shape_}, reverse_shape{shape_}
    { }

    void compute(text_t const & text)
    {
        uint64_t text_length = seqan2::length(text);

        forward_hashes.clear();
        reverse_hashes.clear();
        minimiser_hash.clear();
        minimiser_begin.clear();
        minimiser_end.clear();

        // Return empty vector if text is shorter than k.
        if (k > text_length)
            return;

        reverse_complement_t rc_text{text}; // This does not copy.
        uint64_t possible_minimisers = text_length > w ? text_length - w + 1u : 1u;
        uint64_t possible_kmers = text_length - k + 1;
        uint64_t kmers_per_window = w - k + 1u;

        // Compute all k-mer hashes for both forward and reverse strand.

        // Helper lambda for xor'ing values.
        auto hash_impl = [this] (uint64_t const val)
        {
            return val ^ seed;
        };

        forward_hashes.reserve(possible_kmers);
        reverse_hashes.reserve(possible_kmers);
        auto it1 = seqan2::begin(text);
        auto it2 = seqan2::begin(rc_text);
        seqan2::hashInit(forward_shape, it1);
        seqan2::hashInit(reverse_shape, it2);

        for (uint64_t i = 0; i < possible_kmers; ++i, ++it1, ++it2)
        {
            forward_hashes.push_back(hash_impl(seqan2::hashNext(forward_shape, it1)));
            reverse_hashes.push_back(hash_impl(seqan2::hashNext(reverse_shape, it2)));
        }

        // Choose the minimisers.
        minimiser_hash.reserve(possible_minimisers);
        minimiser_begin.reserve(possible_minimisers);
        minimiser_end.reserve(possible_minimisers);

        // Stores hash, begin and end for all k-mers in the window
        std::deque<std::tuple<uint64_t, uint64_t, uint64_t>> window_values;

        // Initialisation. We need to compute all hashes for the first window.
        for (uint64_t i = 0; i < kmers_per_window; ++i)
        {
            // Get smallest canonical k-mer.
            uint64_t forward_hash = forward_hashes[i];
            uint64_t reverse_hash = reverse_hashes[possible_kmers - i - 1];
            window_values.emplace_back(std::min(forward_hash, reverse_hash), i, i + k - 1);
        }

        auto min = std::min_element(std::begin(window_values), std::end(window_values));
        minimiser_hash.push_back(std::get<0>(*min));
        minimiser_begin.push_back(std::get<1>(*min));
        minimiser_end.push_back(std::get<2>(*min));

        // For the following windows, we remove the first window k-mer (is now not in window) and add the new k-mer
        // that results from the window shifting
        bool minimiser_changed{false};
        for (uint64_t i = 1; i < possible_minimisers; ++i)
        {
            // Shift the window.
            // If current minimiser leaves the window, we need to decide on a new one.
            if (min == std::begin(window_values))
            {
                window_values.pop_front();
                min = std::min_element(std::begin(window_values), std::end(window_values));
                minimiser_changed = true;
            }
            else
            {
                window_values.pop_front();
            }

            uint64_t forward_hash = forward_hashes[kmers_per_window - 1 + i];
            uint64_t reverse_hash = reverse_hashes[possible_kmers - kmers_per_window - i];
            window_values.emplace_back(std::min(forward_hash, reverse_hash),
                                       kmers_per_window + i - 1,
                                       kmers_per_window + i + k - 2);

            if (std::get<0>(window_values.back()) < std::get<0>(*min))
            {
                min = std::prev(std::end(window_values));
                minimiser_changed = true;
            }

            if (minimiser_changed)
            {
                minimiser_hash.push_back(std::get<0>(*min));
                minimiser_begin.push_back(std::get<1>(*min));
                minimiser_end.push_back(std::get<2>(*min));
                minimiser_changed = false;
            }
        }
        return;
    }
};
#endif // SEQAN3_HAS_SEQAN2
