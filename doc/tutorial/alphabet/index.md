# Alphabets in SeqAn {#tutorial_alphabets}

***Learning Objective:***

In this tutorial we look at alphabets and you will learn how to work with nucleotides and amino acids in SeqAn.
We guide you through the most important properties of SeqAn's alphabets and show you the different implemented types.
After completion, you will be able to use the alphabets inside of STL containers and to compare alphabet values.

\tutorial_head{Easy, 45 min, \ref setup, None}

The links on this page mostly point straight into the API documentation which you should use as a reference.
The code examples and assignments are designed to provide some practical experience with our interface
as well as a code basis for your own program development.

[TOC]

# Introduction

An alphabet is the set of symbols of which a biological sequence – or in general a text – is composed.
SeqAn implements specific and optimised alphabets not only for sequences of RNA, DNA and protein components,
but also for quality, secondary structure and gap annotation as well as combinations of the aforementioned.

\assignment{Task}
Read the section *Detailed Description* of the API reference page for \ref alphabet.
This is a detailed introduction to the alphabet module and demonstrates its main advantages.
\endassignment

# The nucleotide alphabets

Nucleotides are the components of (Deoxy)Ribonucleic acid (DNA/RNA) and contain one of the nucleobases
Adenine (A), Cytosine (C), Guanine (G), Thymine (T, only DNA) and Uracil (U, only RNA).
In SeqAn the alphabets seqan3::dna4 and seqan3::rna4 contain exactly the four respective nucleotides.
The trailed number in the alphabets' name represents the number of entities the alphabet holds –
we denote this number as *alphabet size*.
For instance, the alphabet seqan3::dna5 represents five entities as it contains the additional symbol 'N'
to refer to an unknown nucleotide.

## Construction and assignment of alphabet symbols
Let's look at some example code which demonstrates how objects of the seqan3::dna4 alphabet are assigned
from characters.

\snippet alphabet_main.cpp create
\snippet alphabet_main.cpp closing

We have shown three solutions for assigning variables of alphabet type.
1. Construction by character literal, i.e. appending the operator `_dna4` to the respective char symbol. <br>
   This is the handiest way as it can be also used as a temporary object.
2. Assignment by `char` via the global function seqan3::assign_char. <br>
   This is useful if the assignment target already exists, e.g. in a sequence vector.
3. Assignment by rank via the global function seqan3::assign_rank. <br>
   May be used when the *rank* is known.

## The rank of an alphabet symbol
The rank of a symbol is a number in range \[0..alphabet_size) where each number is paired with
an alphabet symbol by a bijective function. In SeqAn the rank is always determined by the lexicographical
order of the underlying characters. For instance, in seqan3::dna4 the bijection is <br>
<code>'A'_dna4 ⟼ 0</code> <br>
<code>'C'_dna4 ⟼ 1</code> <br>
<code>'G'_dna4 ⟼ 2</code> <br>
<code>'T'_dna4 ⟼ 3</code>.

SeqAn provides the function seqan3::to_rank for converting a symbol to its rank value
as demonstrated in the following code example. Note that the data type of the rank is usually the smallest possible
unsigned type that is required for storing the values of the alphabet.

\snippet alphabet_main.cpp rank

## The char representation of an alphabet symbol

Our alphabets also have a character representation because it is more intuitive to work
with them than using the rank. Each alphabet symbol is represented by its respective character
whenever possible (<code>A ⟼ 'A'</code>). Analogously to the rank, SeqAn provides the function
seqan3::to_char for converting a symbol to its character representation.

\snippet alphabet_main.cpp char

Above you have seen that you can assign an alphabet symbol from a character with seqan3::from_char.
In contrast to the rank interface, this assignment is not a bijection because the whole spectrum of available
chars is mapped to values inside the alphabet. For instance, assigning to seqan3::dna4 from any character other
than `C`, `G` or `T` results in the value <code>'A'_dna4</code> and assigning from any character except
`A`, `C`, `G` or `T` to seqan3::dna5 results in the value <code>'N'_dna5</code>. You can avoid the implicit
conversion by using seqan3::assign_char_strict which throws seqan3::invalid_char_assignment on invalid characters.

\snippet alphabet_main.cpp char_strict

You can test the validity of a character by calling seqan3::char_is_valid_for.
It returns true if the character is valid and false otherwise.

## Obtaining the alphabet size
You can retrieve the alphabet size by accessing the class member variable `alphabet_size`
which is implemented in most seqan3::alphabet instances.

\snippet alphabet_main.cpp size

## containers over alphabets

In SeqAn you can use the STL containers to model e.g. sequences, sets or mappings with our alphabets.
The following example shows some exemplary contexts for their use.
For **sequences** we recommend the std::vector with one of SeqAn's alphabet types.
Please note how easily a sequence can be created via the string literal.

\snippet alphabet_main.cpp containers

## Comparability

All alphabets in SeqAn are regular and comparable. This means that you can
use the `<`, `<=`, `>` and `>=` operators to compare the values based on the rank.
For instance, this is useful for sorting a text over the alphabet. Regular implies that the
equality and inequality of two values can be tested with `==` and `!=`.

\snippet alphabet_main.cpp compare

## Example

To wrap up this section on the nucleotide alphabet, the following assignment lets you
practise the use of a SeqAn alphabet and its related functions.
It will also show you a handy advantage of using a vector over an alphabet instead
of using `std::string`: The rank representation can be used straight as an array
index (opposed to using a map with logarithmic access times, for example).

\assignment{Assignment 1: GC content of a sequence}
An important property of DNA and RNA sequences is the *GC content*,
which is the percentage of nucleobases that are either Guanine or Cytosine.
Given the nucleotide counts \f$n_A\f$, \f$n_T\f$, \f$n_G\f$, \f$n_C\f$ the GC content \f$c\f$ is calculated as
\f[ c = \frac{n_G + n_C}{n_A + n_T + n_G + n_C} \f]
Write a program that
1. reads a sequence as command line argument into a vector of seqan3::dna5,
2. counts the number of occurrences for each nucleotide in an array of size `alphabet size` and
3. calculates the GC content.

The seqan3::dna5 type ensures that invalid characters in the input sequence are converted to
<code>'N'</code>. Note that these characters should not influence the GC content.

Pass the sequences `CATTACAG` (3/8 = 37.5%) and `ANNAGAT` (1/5 = 20%) to your program and check
if your results are correct.

\endassignment
\solution
\snippet alphabet_gc_content.cpp exercise
\endsolution

# Other alphabets

Until now, we have focused on alphabets for nucleotides to introduce the properties of SeqAn's alphabet
on a specific example. SeqAn implements, however, many more alphabets. In this section, we want to give you
an overview of the other existing alphabets.

## The amino acid alphabet

Proteins consist of one or more long chains of amino acids. The so-called primary structure of a protein is
expressed as sequences over an amino acid alphabet. The seqan3::aa27 alphabet contains the standard one-letter code
of the 20 canonical amino acids as well as the two proteinogenic amino acids, a termination symbol and
some wildcard characters. For details read the \ref aminoacid page.

## Structure and quality alphabets

The alphabets for structure and quality are sequence *annotations* since they describe additional
properties of the respective sequence.
We distinguish between three types:
1. **Quality alphabet for nucleotides**. The values are produced by sequencing machines and represent the probability
   that a nucleobase was recorded incorrectly. The characters are most commonly found in FASTQ files.
   See \ref quality for details.
2. **RNA structure alphabets**. They describe RNA nucleobases as unpaired or up-/downstream paired and can be found
   in annotated RNA sequence and alignment files (e.g. Stockholm format). Currently we provide the
   [Dot Bracket](\ref seqan3::dot_bracket3) and [WUSS](\ref seqan3::wuss) formats.
3. **Protein structure alphabet**. The [DSSP](\ref seqan3::dssp9) format represents secondary structure elements like
   alpha helices and turns.

You can build an [Alphabet Tuple Composite](\ref seqan3::alphabet_tuple_base) with a nucleotide and quality
alphabet, or nucleotide / amino acid and structure alphabet that stores both information together.
For the use cases just described we offer pre-defined composites (seqan3::qualified, seqan3::structured_rna,
seqan3::structured_aa). See our API documentation for a detailed description of each.

## Gap alphabet

The seqan3::gap alphabet is the smallest alphabet in SeqAn, consisting only of the gap character.
It is most often used in an [Alphabet Variant](\ref seqan3::alphabet_variant) with a nucleotide or amino acid alphabet
to represent gapped sequences, e.g. in alignments. To create a gapped alphabet simply use seqan3::gapped<> with
the alphabet type you want to refine.

\snippet alphabet_main.cpp gapped
