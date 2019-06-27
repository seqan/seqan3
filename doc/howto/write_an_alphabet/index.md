# How to write your own alphabet {#howto_write_an_alphabet}

[TOC]

This HowTo documents how to write a custom alphabet that can be used with the algorithms and data structures in SeqAn3.

\tutorial_head{Moderate, 45 min, \ref tutorial_concepts\, \ref tutorial_alphabets, \ref alphabet}

# Motivation

In the [alphabet tutorial](\ref tutorial_alphabets) you have learned that alphabets are a core component in SeqAn
and represent the smallest unit of biological sequence data. We introduced the common alphabets for
nucleotide, amino acid and gap as well as structure and quality annotation.
However, in your SeqAn application you may want to implement a custom alphabet in order to work efficiently
with SeqAn's algorithms. To achieve this, your custom alphabet must meet certain requirements, which are defined
in [concepts](\ref tutorial_concepts).

For detailed information on the alphabet concepts please read the \ref alphabet page.
In the following sections we demonstrate how to write an alphabet that models them.

A brief summary of the concepts used in this HowTo:
- seqan3::Semialphabet requires your type to have a numerical representation (rank), as well as
  an alphabet size member and comparison operators.
- seqan3::WritableSemialphabet additionally allows to assign a numerical value to an object of your Semialphabet.
- seqan3::Alphabet requires your type to produce a visual representation, which is usually a certain character
  that can be printed as part of a sequence. An Alphabet must also model seqan3::Semialphabet.
- seqan3::WritableAlphabet additionally allows to assign the visual representation to an object of your Alphabet.

# Step by step: Create your own alphabet

In the alphabet tutorial we have calculated the GC content of a nucleotide sequence.
Guanine and Cytosine are complementary nucleobases,
which pair in a DNA molecule by building 3 hydrogen bonds. Adenine and Thymine pair with only 2 hydrogen bonds.
As a consequence, we denote Guanine and Cytosine as strong (S) and Adenine and Thymine as weak (W) nucleobases.
In this section we want to implement a seqan3::Alphabet that consists of the characters `S` and `W` to represent
strong and weak nucleobases.

Let's start with a simple struct that only holds the alphabet's numerical representation, namely the **rank** value:
\snippet dna2_only_rank.cpp struct

If you want SeqAn's algorithms to accept it as an alphabet, you need to make sure that your type
satisfies the requirements of the seqan3::Alphabet concept. A quick check can reveal that this is not the case:
\snippet dna2_only_rank.cpp alphabet_concept

\note
A `static_assert()` will only tell you whether the expression is true or not.
If you want the compiler to tell you **why** the concept is not modelled,
you can construct a dummy requirement like this:
\snippet dna2_alphabet.cpp dummy_requirement
This will fail with a (slightly long and need-to-get-used-to) error message telling you
that `foo()` cannot be called with `dna2` because `constraints are not satisfied ... `.

## Prerequisites

A look at the documentation of seqan3::Alphabet will reveal that it is actually a refinement of other concepts,
more precisely seqan3::Semialphabet which in turn refines std::CopyConstructible and std::StrictTotallyOrdered.
Let's check those:
\snippet dna2_only_rank.cpp other_concepts

You should see that your type models only the std::CopyConstructible concept.
Let's have a look at the documentation for std::StrictTotallyOrdered. Can you guess what it describes?

It describes the requirements for types that are comparable via `<`, `<=`, `>` and `>=`.
This is useful so that you can sort a text over the alphabet, for example.
Additionally, std::StrictTotallyOrdered is again a refinement of std::EqualityComparable,
which requires the `==` and `!=` operators, so let's first implement those with free member functions.
\snippet dna2_semialphabet.cpp equality

\assignment{Excercise}
Implement the inequality operator (!=) for `dna2`.
\endassignment
\solution
\snippet dna2_semialphabet.cpp inequality
\endsolution

We see that our type now models std::EqualityComparable, which is a prerequisite of std::StrictTotallyOrdered.
\snippet dna2_comparison_operators.cpp equality_comparable_concept

\assignment{Excercise}
Implement the four comparison operators *as free member functions* and verify that your type models
std::StrictTotallyOrdered.
\endassignment
\solution
\snippet dna2_comparison_operators.cpp comparison
\endsolution

## Semialphabet

Let's move on to the more interesting concepts. seqan3::Semialphabet requires the *rank interface*
that we have introduced in the [alphabet tutorial](\ref tutorial_alphabets). If we define the member functions
`to_rank`, `assign_rank` and `alphabet_size` for our type, SeqAn automatically exposes them to the respective global
functions `seqan3::to_rank`, `seqan3::assign_rank_to` and `seqan3::alphabet_size`.
Note that `assign_rank` is only required for seqan3::WritableSemialphabet.
\snippet dna2_semialphabet.cpp semialphabet

As you can see from the `static_assert` our dna2 alphabet now models seqan3::Semialphabet and
seqan3::WritableSemialphabet:
\snippet dna2_semialphabet.cpp writable_semialphabet_concept

With this we can automatically use some free functions on our dna2 alphabet. If you are interested in how this works,
read our section on [Customisation points](\ref about_customisation).

## Alphabet

Now that you have a feeling for concepts, have a look at seqan3::Alphabet and seqan3::WritableAlphabet and make
your type also model these concepts. If you keep the suggested namings, SeqAn will automatically
expose the members to the respective global functions that model the concept.

\assignment{Excercise}
1. Implement the member function `char to_char(void)` which returns the char 'S' or 'W' depending on the rank value.
2. Implement the member function `dna2 & assign_char(char const)`.
3. seqan3::WritableAlphabet additionally requires that the type has a member function `bool char_is_valid`
   which can be used to determine whether a given character is an actual representation of the alphabet.
   For the `dna2` alphabet 'S' and 'W' as well as 's' and 'w' should be valid, but nothing else.
   <br>
   Implement the static member function `bool char_is_valid(char const)`.
\endassignment
\solution
\snippet dna2_alphabet.cpp writable_alphabet
\endsolution

At this point the seqan3::Alphabet concept should be modelled successfully and even seqan3::WritableAlphabet 
is fine because we implemented `assign_char` and `char_is_valid`.
\snippet dna2_alphabet.cpp writable_alphabet_concept

## Shortcut: Alphabet base template

In reality, you do not need to define all the functions you have learned in this exercise manually.
Instead, you can inherit your type from seqan3::alphabet_base and just define the char-rank conversion
in both directions. Read the documentation of seqan3::alphabet_base for details and examples.

# Further examples
## Implementation as enum class

This is an example of a minimal custom alphabet that provides implementations for all necessary customisation
points. As an enum class it already models std::StrictTotallyOrdered.
\snippet test/unit/alphabet/custom_alphabet_test.cpp my_alph

## Implementation of a non-default-constructible class

This is an example of a custom alphabet that is not default-constructible and that has a non-default overload for
seqan3::char_is_valid_for.
\snippet test/unit/alphabet/custom_alphabet2_test.cpp my_alph

\note
You should really make your alphabet types [no-throw-default-constructible](\ref std::is_nothrow_default_constructible)
if you can!
