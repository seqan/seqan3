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
- seqan3::Semialphabet requires your type to produce a numerical representation.
  A Semialphabet has an alphabet size and comparison operators as well.
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

Let's start with a simple struct that only holds the alphabet's rank value: 
~~~cpp
#include <seqan3/alphabet/concept.hpp>          // for seqan3::Alphabet concept checks

struct dna2
{
    uint8_t rank{};
};
~~~

If you want SeqAn's algorithms to accept it as an alphabet, you need to make sure that your type
satisfies the requirements of the seqan3::Alphabet concept. A quick check can reveal that this is not the case:
~~~cpp
static_assert(seqan3::Alphabet<dna2>);          // compiler error
~~~

## Prerequisites

A look at the documentation of seqan3::Alphabet will reveal that it is actually a refinement of other concepts, 
more precisely seqan3::Semialphabet which in turn refines std::CopyConstructible and std::StrictTotallyOrdered. 
Let's check those:

~~~cpp
static_assert(std::CopyConstructible<dna2>);    // OK
static_assert(std::StrictTotallyOrdered<dna2>); // compiler error
static_assert(seqan3::Semialphabet<dna2>);      // compiler error
static_assert(seqan3::Alphabet<dna2>);          // compiler error
~~~

You should see that your type models only the std::CopyConstructible concept.
Let's have a look at the documentation for std::StrictTotallyOrdered. Can you guess what it describes?

It describes the requirements for types that are comparable via `<`, `<=`, `>` and `>=`.
This is useful so that you can sort a text over the alphabet, for example.
Additionally, in order to model std::StrictTotallyOrdered we need to make our alphabet
[equality comparable](\ref std::EqualityComparable) by implementing the (in-)equality operators (as free functions).
\snippet alphabet_custom_steps.cpp equality

\assignment{Excercise}
Implement the inequality operator (!=) for `dna2`.
\endassignment
\solution
\snippet alphabet_custom_steps.cpp inequality
\endsolution

We see that our type now models std::EqualityComparable, which is a prerequisite of std::StrictTotallyOrdered.
\snippet alphabet_custom_steps.cpp equality_comparable

\assignment{Excercise}
Implement the four compare operators *as free functions* and verify that your type models std::StrictTotallyOrdered.
\endassignment
\solution
\snippet alphabet_custom_steps.cpp comparison
\endsolution

## Semialphabet

Let's move on to the more interesting concepts. seqan3::Semialphabet requires the *rank* interface
that we have introduced in the [alphabet tutorial](\ref tutorial_alphabets). If we define the member functions
`to_rank`, `assign_rank` and `alphabet_size` for our type, SeqAn automatically exposes them to the respective global
functions `seqan3::to_rank`, `seqan3::assign_rank_to` and `seqan3::alphabet_size`. 
Note that `assign_rank` is only required for seqan3::WritableSemialphabet.
\snippet alphabet_custom_steps.cpp writable_semialphabet

## Alphabet

Now that you have a feeling for concepts, have a look at seqan3::Alphabet and make your type 
also model the seqan3::Alphabet concept. If you keep the suggested namings, SeqAn will automatically expose
the members to the respective global functions that model the concept.

\assignment{Excercise}
Please implement the member function `to_char(void)` which returns the char 'S' or 'W'
depending on the rank value.
\endassignment
\solution
\snippet alphabet_custom_steps.cpp alphabet
\endsolution

\assignment{Excercise}
Now model the full seqan3::WritableAlphabet concept.

1. Implement the member function `dna2 & assign_char(char const)`.
2. Implement the static member function `bool char_is_valid(char const)`.
\endassignment
\solution
\snippet alphabet_custom_steps.cpp writable_alphabet
\endsolution

At this point the seqan3::WritableAlphabet concept should be modelled successfully:
\snippet alphabet_custom_steps.cpp writable_alphabet_concept

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
