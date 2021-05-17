// generated from doc/fragments/@target_alphabet@_implicit_conversion_from_@source_alphabet@.hpp.in

namespace seqan3::doxygen
{
/*!\details
 *
 * Normally, we do not allow implicit conversion of single argument constructors, but in this case we make an exception,
 * because seqan3::rna5 and seqan3::dna5 are interchangeable as they behave nearly the same (e.g. same ranks, same
 * char to rank conversion).
 *
 * \snippet test/snippet/alphabet/nucleotide/rna5_implicit_conversion_from_dna5.cpp main
 *
 * `seqan3::sequence`s (e.g. seqan3::rna5_vector) in general aren't implicitly convertible and must be explicitly
 * copied to be converted:
 *
 * \snippet test/snippet/alphabet/nucleotide/rna5_implicit_conversion_from_dna5_vector.cpp main
 *
 * You can avoid this copy by using `std::ranges::view`s:
 *
 * \snippet test/snippet/alphabet/nucleotide/rna5_implicit_conversion_from_dna5_views.cpp main
 *
 * This conversion constructor only allows converting seqan3::dna5 to seqan3::rna5. Other alphabets that inherit
 * from seqan3::dna5 will not be implicitly convertible to seqan3::rna5.
 *
 * \snippet test/snippet/alphabet/nucleotide/rna5_implicit_conversion_from_dna5_inherit.cpp main
 */
using rna5_implicit_conversion_from_dna5 = void;
}
