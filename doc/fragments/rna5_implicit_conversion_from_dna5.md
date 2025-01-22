<!-- SPDX-FileCopyrightText: 2006-2025 Knut Reinert & Freie Universität Berlin
     SPDX-FileCopyrightText: 2016-2025 Knut Reinert & MPI für molekulare Genetik
     SPDX-License-Identifier: CC-BY-4.0
-->

Normally, we do not allow implicit conversion of single argument constructors, but in this case we make an exception,
because seqan3::rna5 and seqan3::dna5 are interchangeable as they behave nearly the same (e.g. same ranks, same
char to rank conversion).
<br>
\snippet test/snippet/alphabet/nucleotide/rna5_implicit_conversion_from_dna5.cpp main

<br>
`seqan3::sequence`s (e.g. seqan3::rna5_vector) in general are not implicitly convertible and must be explicitly
copied to be converted:
<br>
\snippet test/snippet/alphabet/nucleotide/rna5_implicit_conversion_from_dna5_vector.cpp main

<br>
You can avoid this copy by using `std::ranges::view`s:
<br>
\snippet test/snippet/alphabet/nucleotide/rna5_implicit_conversion_from_dna5_views.cpp main

<br>
This conversion constructor only allows converting seqan3::dna5 to seqan3::rna5. Other alphabets that inherit
from seqan3::dna5 will not be implicitly convertible to seqan3::rna5.
<br>
\snippet test/snippet/alphabet/nucleotide/rna5_implicit_conversion_from_dna5_inherit.cpp main
