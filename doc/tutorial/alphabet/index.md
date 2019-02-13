# Alphabets in SeqAn3 {#tutorial_alphabets}

***Learning Objective:***

In this tutorial we look at alphabets and you will learn how to work with nucleotides and amino acids in SeqAn3.
We guide you through the most important properties of SeqAn's alphabets and show you related concepts.
After completion you will be able to use the alphabets inside of STL containers and to implement your own alphabet.

\tutorial_head{Moderate, 60 min, \ref setup, [Concepts](https://en.cppreference.com/w/cpp/language/constraints)}

The links on this page mostly point straight into the API documentation, which you should use as a reference.
The code examples and exercises are designed to provide some practical experience with our interface
as well as code basis for your own program development.

[TOC]

<br/>

---

# Introduction

An alphabet is the set of symbols, of which a biological sequence – or in general a text – is composed.
SeqAn implements specific and optimised alphabets not only for sequences of RNA, DNA and protein components,
but also for quality, secondary structure and gap annotation, as well as combinations of the aforementioned.

\assignment{Task}
Read the section *Detailed Description* of the API reference page for \ref alphabet.
This is a detailed introduction to the Alphabet module and demonstrates its main advantages.
\endassignment

<br/>

---

# The nucleotide alphabets

Nucleotides are the components of (Deoxy)Ribonucleic acid (DNA/RNA) and contain one of the nucleobases
Adenine (A), Cytosine (C), Guanine (G), Thymine (T, only DNA) and Uracil (U, only RNA).
In SeqAn the alphabets seqan3::dna4 and seqan3::rna4 contain exactly the four respective nucleotides.
The trailed number in the alphabets' name represents the number of entities the alphabet holds –
we denote this number as *alphabet size*.
For instance, the alphabet seqan3::dna5 contains the additional symbol 'N' to refer to an unknown nucleotide.

<br/>

## Creation of alphabet symbols
Let's look at an example code, which demonstrates how characters of the seqan3::dna4 alphabet are created.
\snippet alphabet_main.cpp create
\snippet alphabet_main.cpp closing

We have shown three solutions for assigning variables of alphabet type.
1. Assignment by character literal, i.e. appending the operator `_dna4` to the respective char symbol. <br/>
   This is the handiest way, as it can be used temporary (without creating a variable) as well.
2. Assignment by char via the global function seqan3::assign_char. <br/>
   This is useful if the assignment target already exists, e.g. in a sequence vector.
3. Assignment by rank via the global function seqan3::assign_rank. <br/>
   May be used when the *rank* is known.

<br/>

## The rank of a symbol
The rank of a symbol is a number in range \[0..alphabet_size), where each number is paired with
an alphabet symbol by a bijective function. For instance in seqan3::dna5 the bijection is
<br/> `A ⟼ 0` <br/> `C ⟼ 1` <br/> `G ⟼ 2` <br/> `T ⟼ 3` <br/> `N ⟼ 4`.

SeqAn provides the function seqan3::to_rank for converting a symbol to its rank value,
as demonstrated in the following code example. Note that the data type of the rank is usually the smallest possible
unsigned type that is required for storing the values of the alphabet.
\snippet alphabet_main.cpp rank

<br/>

## Obtaining the alphabet size
You can retrieve the alphabet size by accessing the class member variable `value_size`,
which is implemented (but not required) in all seqan3::Alphabet instances.
However, through the concept it is guaranteed that a global seqan3::alphabet_size_v constant exists.
\snippet alphabet_main.cpp size

<br/>

---

# The Alphabet and NucleotideAlphabet concepts

Concepts are compile-time constraints for template parameters. They are boolean predicates that limit the
set of accepted template types or values and can therefore specialise a function or class. Let's look at an example:

~~~cpp
template <typename alphabet_type>
alphabet_type complement(alphabet_type const alph)
{
    // implementation...
}
~~~

We want to implement a function `complement` that computes the complementary nucleobase of a given one.
Because we want a unique interface for each alphabet, e.g. seqan3::dna4 and seqan3::dna5, we implement
it as a template function as shown above. The problem is that in this way the function does not only accept DNA
and RNA alphabets – for which this function is useful – but also any other type, like `float` or `std::string`.

The solution is to *constrain* the function to a certain set of alphabets, namely nucleotide alphabets.
To achieve this, SeqAn defines a concept seqan3::NucleotideAlphabet, where is specified how nucleotide
alphabets look like, so the compiler can decide, whether a given type is valid or not. In our example
we can now require that the given `alphabet_type` is a nucleotide alphabet.

~~~cpp
template <typename alphabet_type>
    requires seqan3::NucleotideAlphabet<alphabet_type>
alphabet_type complement(alphabet_type const alph)
{
    // implementation...
}
~~~

If you try to call the function `complement` on e.g. a `float` variable, the compilation fails with a
helpful error message stating that `float` is not a seqan3::NucleotideAlphabet.

A more general concept in SeqAn is seqan3::Alphabet. It ensures for all alphabet types in SeqAn that you can
- get and assign the rank (seqan3::to_rank, seqan3::assign_rank)
- get and assign the char representation (seqan3::to_char, seqan3::assign_char, seqan3::assign_char_strict)
- check if a char is valid in this alphabet (seqan3::char_is_valid_for)
- retrieve the alphabet size (seqan3::Alphabet::alphabet_size_v) and
- compare their values (operators `==`, `>`, `>=`, ...).

\note Of course, the seqan3::NucleotideAlphabet satisfies seqan3::Alphabet. <br/>
      It actually extends the alphabet concept with the [complement()](\ref seqan3::NucleotideAlphabet::complement)
      function.

<br/>

---

# Containers over alphabets

In SeqAn you can use the STL containers to model e.g. sequences, sets or mappings with our alphabets.
For sequences we recommend the std::vector with one of SeqAn's alphabet types. Instead of building a map from
characters with logarithmic access times, the rank representation can be used straight as an array index
with constant access time. In the following exercise you can practise the use of alphabet containers together
with the previously discussed functions.

\assignment{Excercise: GC content of a sequence}
An important property of DNA and RNA molecules is the *GC content*,
which is the percentage of nucleobases that are either Guanine or Cytosine.
Given the nucleotide counts \f$n_A\f$, \f$n_T\f$, \f$n_G\f$, \f$n_C\f$ the GC content \f$c\f$ is calculated as
\f[ c = \frac{n_G + n_C}{n_A + n_T + n_G + n_C} \f]
Write a program that
1. reads a sequence as command line argument into a vector of seqan3::dna4,
2. fails for incorrect characters in the sequence,
3. counts the number of occurrences for each nucleotide and
4. calculates the GC content.
\endassignment
\solution
\snippet alphabet_gc_content.cpp exercise
\endsolution

<br/>

---

# Other alphabets

Until now, we have focused on alphabets for nucleotides to introduce the properties of SeqAn's alphabet
on a specific example. SeqAn implements, however, many more alphabets. In this section, we want to give you
an overview to the existing alphabets and in the end you can implement an additional one yourself.

<br/>

## The amino acid alphabet

Proteins consist of one or more long chains of amino acids. The so-called primary structure of a protein is
expressed as sequences over an amino acid alphabet. The seqan3::aa27 alphabet contains the standard one-letter code
of the 20 canonical amino acids, as well as the two proteinogenic amino acids, a termination symbol and
some wildcard characters. For details read the seqan3::Aminoacid page.

<br/>

## Structure and quality alphabets

The alphabets for structure and quality are usually sequence *annotations*, as they describe additional
properties of the respective sequence. Therefore, they usually build a
[Cartesian Composition](\ref seqan3::cartesian_composition) with a nucleotide or amino acid alphabet.
We distinguish three types:
1. Quality alphabet for nucleotides. The values are produced by sequencing machines and represent the probability
   that a nucleobase was recorded incorrectly. The characters are most commonly found in FASTQ files.
   See seqan3::Quality for details.
2. RNA structure alphabets. They describe RNA nucleobases as unpaired or up-/downstream paired and can be found
   in annotated RNA sequence and alignment files (e.g. Stockholm format). Currently we provide the
   [Dot Bracket](\ref seqan3::dot_bracket3) and [WUSS](\ref seqan3::wuss) formats.
3. Protein structure alphabet. The [DSSP](\ref seqan3::dssp9) format represents secondary structure elements like
   alpha helices and turns.

<br/>

## Gap alphabet

The seqan3::gap alphabet is the smallest alphabet in SeqAn, consisting of the gap character only.
It is used in a [Union Composition](\ref seqan3::union_composition) with a nucleotide or amino acid alphabet
to represent gapped sequences, e.g. in alignments. To create a gapped alphabet simply use seqan3::gapped<> with
the alphabet type you want to extend.
\snippet alphabet_main.cpp gapped

<br/>

## Create your own alphabet

Finally, we want to demonstrate how easy it is to implement a custom alphabet with SeqAn3.

In the previous exercise we have calculated the GC content. Guanine and Cytosine are complementary nucleobases,
which pair in a DNA molecule by building 3 hydrogen bonds. Adenine and Thymine pair with only 2 hydrogen bonds.
As a consequence, we denote Guanine and Cytosine as strong (S) and Adenine and Thymine as weak (W) nucleobases.

\assignment{Excercise: A little alphabet}

In this exercise we want to implement a seqan3::Alphabet that consists of the characters `S` and `W` to represent
strong and weak nucleobases. Therefore, we have to implement all requirements for the alphabet concept,
which are documented [here](\ref seqan3::Alphabet). As you see, a concept can also hierarchically require
other concepts, e.g. a seqan3::Alphabet has to satisfy seqan3::Semialphabet as well.

Here is a code frame and a test function to start with:
~~~cpp
#include <iostream>   // std::cerr, std::endl
#include <seqan3/alphabet/all.hpp>

using namespace seqan3;

class dna2
{
    // ... implementation of member functions and variables
}

// ... implementation of operators for std::EqualityComparable
// ... implementation of operators for std::StrictTotallyOrdered

// Constrained function that works only for seqan3::Alphabet types.
template <typename alph_type>
    requires Alphabet<alph_type>
void test_function(alph_type)
{
    std::cerr << "You're good!" << std::endl;
    std::cerr << "The alphabet size is " << (unsigned)alphabet_size_v<dna2> << "." << std::endl;
}

int main ()
{
    // Let's test our new alphabet class here. The compilation fails, if members are missing.
    test_function(dna2{});
    return 0;
}
~~~

The following steps guide you through the implementation.
If you keep the suggested naming, SeqAn will automatically expose the members to the respective global functions.
1. Define the types `rank_type` and `char_type` for your alphabet.
2. Define the static member constant `value_size`, which holds the alphabet size.
3. Define a (private) member to store the rank value.
4. Implement the member functions `to_rank(void)`, `to_char(void)`, `assign_rank(rank_type const)` and
   `assign_char(char_type const)`.
5. Implement the static member function `char_is_valid(char_type const)`.
6. Implement the member function `assign_char_strict(char_type const)` that throws seqan3::invalid_char_assignment,
   if `char_is_valid` evaluates to false.
7. Overload all [equality](https://en.cppreference.com/w/cpp/concepts/EqualityComparable) and
   [comparison](https://en.cppreference.com/w/cpp/concepts/StrictTotallyOrdered) operators for your class.
\endassignment
\solution
\snippet alphabet_custom.cpp exercise
\endsolution
