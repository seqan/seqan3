#include <seqan3/io/detail/output_iterator_conversion_adaptor.hpp>
#include <seqan3/alphabet/nucleotide/dna4.hpp>

using namespace seqan3;

int main()
{
//! [usage]
using insert_iter = ranges::back_insert_iterator<std::vector<dna4>>;
using out_iter    = seqan3::detail::output_iterator_conversion_adaptor<insert_iter, dna4>;

std::vector<dna4> out_vec;
out_iter it{insert_iter{out_vec}};

*it = 'A';
*it = 'C';
*it = 'G';
*it = 'T';

for (auto && c : out_vec)
    std::cout << to_char(c);
std::cout << '\n'; // prints ACGT
//! [usage]
}
