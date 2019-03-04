# Alphabets in SeqAn3 {#tutorial_alphabets}

***Learning Objective:***

In this tutorial we look at alphabets and you will learn how to work with nucleotides and amino acids in SeqAn3.
We guide you through the most important properties of SeqAn's alphabets and show you related concepts.
After completion, you will be able to use the alphabets inside of STL containers and to implement your own alphabet.

\tutorial_head{Moderate, 60 min, \ref setup, [Concepts](https://en.cppreference.com/w/cpp/language/constraints)}

The links on this page mostly point straight into the API documentation, which you should use as a reference.
The code examples and exercises are designed to provide some practical experience with our interface
as well as a code basis for your own program development.

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
For instance, the alphabet seqan3::dna5 represents five entities as it contains the additional symbol 'N'
to refer to an unknown nucleotide.

<br/>

## Creation of alphabet symbols
Let's look at some example code, which demonstrates how characters of the seqan3::dna4 alphabet are created.
\snippet alphabet_main.cpp create
\snippet alphabet_main.cpp closing

We have shown three solutions for assigning variables of alphabet type.
1. Assignment by character literal, i.e. appending the operator `_dna4` to the respective char symbol. <br/>
   This is the handiest way, as it can be used temporarily (without creating a variable) as well.
2. Assignment by `char` via the global function seqan3::assign_char. <br/>
   This is useful if the assignment target already exists, e.g. in a sequence vector.
3. Assignment by rank via the global function seqan3::assign_rank. <br/>
   May be used when the *rank* is known.

<br/>

## The rank of an alphabet symbol
The rank of a symbol is a number in range \[0..alphabet_size), where each number is paired with
an alphabet symbol by a bijective function. For instance in seqan3::dna4 the bijection is
<br/> `A ⟼ 0` <br/> `C ⟼ 1` <br/> `G ⟼ 2` <br/> `T ⟼ 3`.

SeqAn provides the function seqan3::to_rank for converting a symbol to its rank value,
as demonstrated in the following code example. Note that the data type of the rank is usually the smallest possible
unsigned type that is required for storing the values of the alphabet.
\snippet alphabet_main.cpp rank

<br/>

## The char representation of an alphabet symbol

Our alphabets also have a character representation, because it is more intuitive to work
with them than using the rank. Each alphabet symbol is represented by its respective character
whenever possible (`A ⟼ 'A'`). Analogously to the rank, SeqAn provides the function 
seqan3::to_char for converting a symbol to its character representation.

\snippet alphabet_main.cpp char

Above you have seen that you can assign an alphabet symbol from a character with seqan3::from_char. 
In contrast to the rank interface, this assignment is not a bijection, because the whole spectrum of available
chars is mapped to values inside the alphabet. For instance, assigning to seqan3::dna4 from any character other
than `C`, `G` or `T` results in the value 'A'_dna4 and assigning from any character except `A`, `C`, `G` or `T`
to seqan3::dna5 results in the value 'N'_dna5. You can avoid the implicit conversion by using
seqan3::assign_char_strict, which throws seqan3::invalid_char_assignment on invalid characters.

\snippet alphabet_main.cpp char_strict

<br/>

## Obtaining the alphabet size
You can retrieve the alphabet size by accessing the class member variable `value_size`,
which is implemented in all seqan3::Alphabet instances.
\snippet alphabet_main.cpp size

<br/>

---

## Containers over alphabets

In SeqAn you can use the STL containers to model e.g. sequences, sets or mappings with our alphabets.
The following example shows some exemplary contexts for their use. 
For **sequences** we recommend the std::vector with one of SeqAn's alphabet types. 
Please note how easily a sequence can be created via the string literal.
\snippet alphabet_main.cpp containers

<br/>

---

## Example

To wrap up this section on the nucleotide alphabet, the following exercise will let you 
practise the use of a SeqAn alphabet and its related functions.
It will also show you a handy advantage of using a vector over an alphabet instead 
of using `std::string`: The rank representation can be used straight as an array
index (opposed to e.g. using a map with logarithmic access times).

\assignment{Excercise: GC content of a sequence}
An important property of DNA and RNA molecules is the *GC content*,
which is the percentage of nucleobases that are either Guanine or Cytosine.
Given the nucleotide counts \f$n_A\f$, \f$n_T\f$, \f$n_G\f$, \f$n_C\f$ the GC content \f$c\f$ is calculated as
\f[ c = \frac{n_G + n_C}{n_A + n_T + n_G + n_C} \f]
Write a program that
1. reads a sequence as command line argument into a vector of seqan3::dna5,
2. counts the number of occurrences for each nucleotide in an array of size `alphabet size` and
3. calculates the GC content.

The seqan3::dna5 type ensures that invalid characters in the input sequence are converted to 'N'.
Note that these characters should not influence the GC content.
\endassignment
\solution
\snippet alphabet_gc_content.cpp exercise
\endsolution

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
To achieve this, SeqAn defines a concept seqan3::NucleotideAlphabet that specifies how nucleotide
alphabets look like, so the compiler can decide, whether a given type is valid or not. In our example
we can now require that the given `alphabet_type` is a nucleotide alphabet.

~~~cpp
template <seqan3::NucleotideAlphabet alphabet_type>
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

The alphabets for structure and quality are sequence *annotations*, as they describe additional
properties of the respective sequence. 
We distinguish three types:
1. Quality alphabet for nucleotides. The values are produced by sequencing machines and represent the probability
   that a nucleobase was recorded incorrectly. The characters are most commonly found in FASTQ files.
   See seqan3::Quality for details.
2. RNA structure alphabets. They describe RNA nucleobases as unpaired or up-/downstream paired and can be found
   in annotated RNA sequence and alignment files (e.g. Stockholm format). Currently we provide the
   [Dot Bracket](\ref seqan3::dot_bracket3) and [WUSS](\ref seqan3::wuss) formats.
3. Protein structure alphabet. The [DSSP](\ref seqan3::dssp9) format represents secondary structure elements like
   alpha helices and turns.

You can build a [Cartesian Compositions](\ref seqan3::cartesian_composition) with a nucleotide and quality 
alphabet, or nucleotide / amino acid and structure alphabet that stores both information together. 
For the use cases just described we offer pre-defined composites: (list qualified, structured_aa, structured_rna). 
See our API documentation for a detailed description of each.

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
In this section we want to implement a seqan3::Alphabet that consists of the characters `S` and `W` to represent
strong and weak nucleobases. 

Let's start with a simple struct that only holds the alphabet's rank value: 
\snippet alphabet_dna2_struct.cpp struct

But if you want SeqAn's algorithms to accept it as an alphabet, you need to make sure that your type
satisfies the requirements of the seqan3::Alphabet concept. A quick check can reveal that this is not the case:
~~~cpp
static_assert(seqan3::Alphabet<dna2>);          // compiler error
~~~

A look at the documentation of seqan3::Alphabet will reveal that it is actually a refinement of other concepts, 
more precisely seqan3::Semialphabet which in turn refines std::Regular and std::StrictTotallyOrdered. 
Let's check those:

~~~cpp
static_assert(std::Regular<dna2>);              // compiler error
static_assert(std::StrictTotallyOrdered<dna2>); // compiler error
static_assert(seqan3::Semialphabet<dna2>);      // compiler error
static_assert(seqan3::Alphabet<dna2>);          // compiler error
~~~

You should see that your type models none of the concepts. Let's change this one by one.
For std::Regular we need to make our alphabet [equality comparable](\ref std::EqualityComparable)
by implementing the (in-)equality operators (as free functions).
\snippet alphabet_dna2_steps.cpp equality

\assignment{Excercise}
Implement the inequality operator (!=) for `dna2`.
\endassignment
\solution
\snippet alphabet_dna2_steps.cpp inequality
\endsolution

We see that our type now models std::Regular.
\snippet alphabet_dna2_steps.cpp regular

Let's have a look at the documentation for std::StrictTotallyOrdered. 
Can you guess what it describes?

It describes the requirements for types that are comparable via `<`, `<=`, `>` and `>=`.
This is useful so that you can sort a text over the alphabet, for example.

\assignment{Excercise}
Implement the four compare operators *as free functions* and verify that your type models std::StrictTotallyOrdered.
\endassignment
\solution
\snippet alphabet_dna2_steps.cpp compare
\endsolution

The next concept is more interesting for us: seqan3::Semialphabet requires the *rank* interface
that we have introduced in the beginning of this tutorial. If we define the member functions `to_rank`,
`assign_rank` and `value_size` for our type, SeqAn automatically exposes them to the respective global
functions `seqan3::to_rank`, `seqan3::assign_rank` and `seqan3::alphabet_size`.
\snippet alphabet_dna2_steps.cpp semialphabet

\assignment{Excercise}
Now that you have a feeling for concepts, have a look at seqan3::Alphabet and make your type 
also model the full seqan3::Alphabet concept. If you keep the suggested naming, SeqAn will automatically expose 
the members to the respective global functions that model the concept.

1. Define the [char_type](\ref seqan3::alphabet_base::char_type) for your alphabet (the type that you want to emit when calling 
   [to_char()](\ref seqan3::alphabet_base::to_char)).
2. Implement the member functions [to_char()](\ref seqan3::alphabet_base::to_char) and 
   [assign_char(char_type const)](\ref seqan3::alphabet_base::assign_char).
3. Implement the static member function [char_is_valid(char_type const)](\ref seqan3::alphabet_base::char_is_valid).
4. Implement the member function [assign_char_strict(char_type const)](\ref seqan3::alphabet_base::assign_char_strict) 
   that throws seqan3::invalid_char_assignment, if [char_is_valid](\ref seqan3::alphabet_base::char_is_valid) 
   evaluates to false.
\endassignment
\solution
\snippet alphabet_dna2_steps.cpp alphabet
\endsolution

Now the seqan3::Alphabet concept should be modelled successfully:
\snippet alphabet_dna2_steps.cpp alphabet_concept

In reality, you do not need to define all the functions you have learned in this exercise manually. 
Instead, you can inherit you type from seqan3::alphabet_base and just define the char-rank conversion
in both directions. Read the documentation of seqan3::alphabet_base fo details. 
