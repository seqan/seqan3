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

//!\brief Whether to use xor for computing the hash value.
enum use_xor : bool
{
    //!\brief Do not use xor.
    no,
    //!\brief Use xor.
    yes
};

template<typename shape_t, use_xor do_xor = use_xor::yes>
struct minimiser
{
private:
    //!\brief The alphabet type.
    using alphabet_t = seqan::Dna;
    //!\brief The text type.
    using text_t = seqan::String<alphabet_t>;
    //!\brief The type of the complemented text.
    using complement_t = seqan::ModifiedString<text_t const, seqan::ModComplementDna>;
    //!\brief The type of the reverse complemented text.
    using reverse_complement_t = seqan::ModifiedString<complement_t, seqan::ModReverse>;
    //!\brief The seqan::Shape type to compute rolling hashes with.
    //using shape_t = seqan::Shape<alphabet_t, seqan::GenericShape>;

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
     * \param[in] w_    The window size.
     * \param[in] k_    The k-mer size.
     * \param[in] seed_ The seed to use. Default: 0x8F3F73B5CF1C9ADE.
     */
    minimiser(window const w_, kmer const k_, uint64_t const seed_ = 0x8F3F73B5CF1C9ADE) :
        w{w_.v}, k{k_.v}, seed{seed_}
    {
        seqan::resize(forward_shape, k);
        seqan::resize(reverse_shape, k);
    }

    /*!\brief Resize the minimiser.
     * \param[in] w_    The new window size.
     * \param[in] k_    The new k-mer size.
     * \param[in] seed_ The new seed to use. Default: 0x8F3F73B5CF1C9ADE.
     */
     void resize(window const w_, kmer const k_, shape_t new_shape, uint64_t const seed_ = 0x8F3F73B5CF1C9ADE)
   {
       w = w_.v;
       k = k_.v;
       seed = seed_;
       forward_shape = new_shape;
       reverse_shape = new_shape;
       //seqan::resize(forward_shape, k);
       //seqan::resize(reverse_shape, k);
   }

    void compute(text_t const & text)
    {
        uint64_t text_length = seqan::length(text);

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

        // Helper lambda for xor'ing values depending on `do_xor`.
        auto hash_impl = [this] (uint64_t const val)
        {
            if constexpr(do_xor)
                return val ^ seed;
            else
                return val;
        };

        forward_hashes.reserve(possible_kmers);
        reverse_hashes.reserve(possible_kmers);
        seqan::hashInit(forward_shape, seqan::begin(text));
        seqan::hashInit(reverse_shape, seqan::begin(rc_text));

        for (uint64_t i = 0; i < possible_kmers; ++i)
        {
            forward_hashes.push_back(hash_impl(seqan::hashNext(forward_shape, seqan::begin(text) + i)));
            reverse_hashes.push_back(hash_impl(seqan::hashNext(reverse_shape, seqan::begin(rc_text) + i)));
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

/*using namespace seqan;

struct Minimiser
{
public:

    // Random, but static value for xor for hashes. Counteracts consecutive minimisers.
    // E.g., without it, the next minimiser after a poly-A region AAAAA would be most likely something like AAAAC.
    uint64_t const seed{0x8F3F73B5CF1C9ADE};
    // Shape for forward hashes
    Shape<Dna, GenericShape> kmerShape;
    // Shape for hashes on reverse complement
    Shape<Dna, GenericShape> revCompShape;
    // k-mer size
    uint16_t k{19};
    // window size
    uint32_t w{25};
    // start positions of minimisers
    std::vector<uint64_t> minBegin;
    // end positions of minimisers
    std::vector<uint64_t> minEnd;

    template<typename TIt>
    inline void hashInit(TIt it)
    {
        seqan::hashInit(kmerShape, it);
    }

    template<typename TIt>
    inline auto hashNext(TIt it)
    {
        return seqan::hashNext(kmerShape, it);
    }

    template<typename TIt>
    inline void revHashInit(TIt it)
    {
        seqan::hashInit(revCompShape, it);
    }

    template<typename TIt>
    inline auto revHashNext(TIt it)
    {
        return seqan::hashNext(revCompShape, it);
    }

    inline auto length()
    {
        return seqan::length(kmerShape);
    }

    inline void resize(uint16_t newKmerSize, uint32_t neww, Shape<Dna, GenericShape> newkmerShape)
    {
        k = newKmerSize;
        w = neww;
        kmerShape = newkmerShape;
        revCompShape = newkmerShape;
    }

    std::vector<uint64_t> getHash(DnaString const & text)
    {
        if (k > seqan::length(text))
            return std::vector<uint64_t> {};

        // Reverse complement without copying/modifying the original string
        typedef ModifiedString<ModifiedString<DnaString const, ModComplementDna>, ModReverse> TRC;
        TRC revComp(text);

        uint64_t possible = seqan::length(text) > w ? seqan::length(text) - w + 1 : 1;
        uint32_t windowKmers = w - k + 1;

        std::vector<uint64_t> kmerHashes;
        // Stores hash, begin and end for all k-mers in the window
        std::deque<std::tuple<uint64_t, uint64_t, uint64_t>> windowValues;
        kmerHashes.reserve(possible);
        minBegin.reserve(possible);
        minEnd.reserve(possible);

        auto it = begin(text);
        auto rcit = begin(revComp);
        hashInit(it);
        revHashInit(rcit);

        // Initialisation. We need to compute all hashes for the first window.
        for (uint32_t i = 0; i < windowKmers; ++i)
        {
            // Get smallest canonical k-mer
            uint64_t kmerHash = hashNext(it) ^ seed;
            uint64_t revcHash = revHashNext(rcit) ^ seed;
            if (kmerHash <= revcHash)
            {
                uint64_t distance = std::distance(begin(text), it);
                windowValues.push_back(std::make_tuple(kmerHash, distance, distance + k - 1));
            }
            else
            {
                uint64_t distance = std::distance(rcit, end(revComp)) - k;
                windowValues.push_back(std::make_tuple(revcHash, distance, distance + k - 1));
            }
            ++it;
            ++rcit;
        }

        auto min = std::min_element(std::begin(windowValues), std::end(windowValues));
        kmerHashes.push_back(std::get<0>(*min));
        minBegin.push_back(std::get<1>(*min));
        minEnd.push_back(std::get<2>(*min));

        // For the following windows, we remove the first window k-mer (is now not in window) and add the new k-mer
        // that results from the window shifting
        for (uint64_t i = 1; i < possible; ++i)
        {
            if (min == std::begin(windowValues))
            {
                windowValues.pop_front();
                min = std::min_element(std::begin(windowValues), std::end(windowValues));
            }
            else
                windowValues.pop_front();

            uint64_t kmerHash = hashNext(it) ^ seed;
            uint64_t revcHash = revHashNext(rcit) ^ seed;
            if (kmerHash <= revcHash)
            {
                uint64_t distance = std::distance(begin(text), it);
                windowValues.push_back(std::make_tuple(kmerHash, distance, distance + k - 1));
            }
            else
            {
                uint64_t distance = std::distance(rcit, end(revComp)) - k;
                windowValues.push_back(std::make_tuple(revcHash, distance, distance + k - 1));
            }
            ++it;
            ++rcit;

            if (std::get<0>(windowValues.back()) < std::get<0>(*min))
                min = std::end(windowValues) - 1;

            kmerHashes.push_back(std::get<0>(*min));
            minBegin.push_back(std::get<1>(*min));
            minEnd.push_back(std::get<2>(*min));
        }

        return kmerHashes;
    }

    std::vector<uint64_t> getMinimiser(Dna5String const & text)
    {
        if (k > seqan::length(text))
            return std::vector<uint64_t> {};

        // Reverse complement without copying/modifying the original string
        typedef ModifiedString<ModifiedString<Dna5String const, ModComplementDna>, ModReverse> TRC;
        TRC revComp(text);

        uint64_t possible = seqan::length(text) > w ? seqan::length(text) - w + 1 : 1;
        uint32_t windowKmers = w - k + 1;

        std::vector<uint64_t> kmerHashes{};
        kmerHashes.reserve(possible); // maybe rather reserve to expected?

        // Stores hash, begin and end for all k-mers in the window
        std::deque<uint64_t> windowValues;

        auto it = begin(text);
        auto rcit = begin(revComp);
        hashInit(it);
        revHashInit(rcit);

        // Initialisation. We need to compute all hashes for the first window.
        for (uint32_t i = 0; i < windowKmers; ++i)
        {
            // Get smallest canonical k-mer
            uint64_t kmerHash = hashNext(it) ^ seed;
            uint64_t revcHash = revHashNext(rcit) ^ seed;
            if (kmerHash <= revcHash)
            {
                windowValues.push_back(kmerHash);
            }
            else
            {
                windowValues.push_back(revcHash);
            }
            ++it;
            ++rcit;
        }

        auto min = std::min_element(std::begin(windowValues), std::end(windowValues));
        kmerHashes.push_back(*min);


        // For the following windows, we remove the first window k-mer (is now not in window) and add the new k-mer
        // that results from the window shifting
        bool minimiser_changed{false};
        for (uint64_t i = 1; i < possible; ++i)
        {
            if (min == std::begin(windowValues))
            {
                windowValues.pop_front();
                min = std::min_element(std::begin(windowValues), std::end(windowValues));
                minimiser_changed = true;
            }
            else
                windowValues.pop_front();

            uint64_t kmerHash = hashNext(it) ^ seed;
            uint64_t revcHash = revHashNext(rcit) ^ seed;
            if (kmerHash <= revcHash)
            {
                windowValues.push_back(kmerHash);
            }
            else
            {
                windowValues.push_back(revcHash);
            }
            ++it;
            ++rcit;

            if (windowValues.back() < *min)
            {
                min = std::end(windowValues) - 1;
                minimiser_changed = true;
            }

            if (minimiser_changed)
            {
                kmerHashes.push_back(*min);
                minimiser_changed = false;
            }
        }

        return kmerHashes;
    }
};*/
#endif // SEQAN3_HAS_SEQAN2
