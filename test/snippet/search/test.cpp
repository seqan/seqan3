#include <seqan3/search/dream_index/dream_index.hpp>
#include <seqan3/io/stream/debug_stream.hpp>
// #include <seqan3/alphabet/nucleotide/dna4.hpp>
// #include <seqan3/range/container/concept.hpp>

// #include <sdsl/int_vector.hpp>
// #include <sdsl/sd_vector.hpp>

using namespace seqan3;

// template <reservable_container_concept T>
// struct S {};

int main()
{
    // S<detail::bitvector<uncompressed>> s;
    // (void) s;
    // // debug_stream << std::boolalpha << container_concept<detail::bitvector<uncompressed>>;
    // /*
    // detail::binning_directory<ibf, uncompressed> bd(64, 512);
    // debug_stream << "Set\n";
    // bd.set(20,5);
    // debug_stream << "Get\n";
    // debug_stream << bd.get(20) << '\n';
    // */
    // // bv[3] = 1;
    // // debug_stream << bv.data << '\n';
    //
    dream_index di(5, 3, 1024);
    di.insert_data(0, "AGCTACGTAAAAAAAAAAAAAAAA"_dna4);
    di.insert_data(1, "AGCTACGAAAAAAAAAAAAAAAAA"_dna4);
    di.insert_data(2, "AGCTACAAAAAAAAAAAAAAAAAA"_dna4);
    debug_stream << di.count("AGCTACGT"_dna4) << '\n';
    debug_stream << di.get_bins("AGCTACGT"_dna4) << '\n';
    debug_stream << di.get_bins("AGCTACGT"_dna4, 100) << '\n';

}
