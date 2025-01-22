# How to write your own alphabet {#howto_write_an_alphabet}

<!-- SPDX-FileCopyrightText: 2006-2025 Knut Reinert & Freie Universität Berlin
     SPDX-FileCopyrightText: 2016-2025 Knut Reinert & MPI für molekulare Genetik
     SPDX-License-Identifier: CC-BY-4.0
-->

[TOC]

This HowTo documents how to write a custom alphabet that can be used with the algorithms and data structures in SeqAn.

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
- seqan3::semialphabet requires your type to have a numerical representation (rank), as well as
  an alphabet size and comparison operators.
- seqan3::writable_semialphabet additionally requires being able to change the value of the object
  via the rank representation.
- seqan3::alphabet requires that your type has a visual representation in addition to the numerical representation.
  Usually this is a character type like char.
- seqan3::writable_alphabet additionally requires being able to change the value of the object
  via the char representation.

# Step by step: Create your own alphabet

In the alphabet tutorial we have calculated the GC content of a nucleotide sequence.
Guanine and Cytosine are complementary nucleobases,
which pair in a DNA molecule by building 3 hydrogen bonds. Adenine and Thymine pair with only 2 hydrogen bonds.
As a consequence, we denote Guanine and Cytosine as strong (S) and Adenine and Thymine as weak (W) nucleobases.
In this section we want to implement a seqan3::alphabet that consists of the characters `S` and `W` to represent
strong and weak nucleobases.

Let's start with a simple struct that only holds the alphabet's numerical representation, namely the **rank** value.
It is not specified how the value of an alphabet is stored internally, however alphabets typically store the rank as a
member variable.
\snippet dna2_only_rank.cpp struct

\note
The type of the member variable is typically the smallest unsigned integer type that can hold the maximum rank.
In our case we have only two values so `bool` would be sufficient.
However `bool` is not actually smaller than `uint8_t` and does not always behave like an unsigned integral type
so we use `uint8_t` here. You only need a different type if your alphabet has an alphabet size >=255.

If you want SeqAn's algorithms to accept it as an alphabet, you need to make sure that your type
satisfies the requirements of the seqan3::alphabet concept. A quick check can reveal that this is not the case:
\snippet dna2_only_rank.cpp alphabet_concept

\note
A `static_assert()` will only tell you whether the expression is true or not.
If you want the compiler to tell you **why** the concept is not modelled,
you can construct a dummy requirement like this:
\snippet dna2_alphabet.cpp dummy_requirement
This will fail with a (slightly long and need-to-get-used-to) error message telling you
that `foo()` cannot be called with `dna2` because `constraints are not satisfied ... `.

## Prerequisites

A look at the documentation of seqan3::alphabet will reveal that it is actually a refinement of other concepts,
more precisely seqan3::semialphabet which in turn refines std::copy_constructible and std::totally_ordered.
Let's check those:
\snippet dna2_only_rank.cpp other_concepts

You should see that your type models only the std::copy_constructible concept.
Let's have a look at the documentation for std::totally_ordered. Can you guess what it describes?

It describes the requirements for types that are comparable via `<`, `<=`, `>` and `>=`.
This is useful so that you can sort a text over the alphabet, for example.
Additionally, std::totally_ordered is again a refinement of std::equality_comparable,
which requires the `==` and `!=` operators, so let's first implement those.
\snippet dna2_equality_operator.cpp equality

As you can see we chose to implement them as friend functions. You can also implement them as free functions in the
namespace of your class, but you should avoid regular member functions, because they result in different conversion
rules for left-hand-side and right-hand-side.

\assignment{Excercise}
Implement the inequality operator (!=) for `dna2`.
\endassignment
\solution
\snippet dna2_inequality_operator.cpp inequality
\endsolution

We see that our type now models std::equality_comparable, which is a prerequisite of std::totally_ordered.

\assignment{Excercise}
Implement the four comparison operators and verify that your type models std::totally_ordered.
\endassignment
\solution
\snippet dna2_comparison_operators.cpp comparison
\endsolution

## semialphabet

Let's move on to the more interesting concepts. seqan3::semialphabet constitutes the *rank interface*
that we introduced in the [alphabet tutorial](\ref tutorial_alphabets). Have a look at the API reference again.
Beyond the conceptional requirements, it also requires that seqan3::alphabet_size and seqan3::to_rank can be
called on your alphabet.

There are different ways to satisfy seqan3::alphabet_size and seqan3::to_rank, have a look at the respective API
reference and also the [documentation on customisation points](\ref about_customisation).

In this case we choose to implement the functionality as member functions:
\snippet dna2_semialphabet.cpp semialphabet

As you can see from the `static_assert` our dna2 alphabet now models seqan3::semialphabet and
seqan3::writable_semialphabet:
\snippet dna2_semialphabet.cpp writable_semialphabet_concept

You can also try out the customisation points directly (they appear like free functions). Here is an example:
\snippet dna2_semialphabet.cpp free_functions

## alphabet

Now that you have a feeling for concepts, have a look at seqan3::alphabet and seqan3::writable_alphabet and make
your type also model these concepts.

\assignment{Excercise}
1. Implement the member function `char %to_char()` which returns the char 'S' or 'W' depending on the rank value.
2. Implement the member function `dna2 & assign_char(char const)`.
3. There is a function object that tells us which characters are valid for a given alphabet: seqan3::char_is_valid_for.
   By default, all characters are "valid" that are preserved when being assigned from and then be converted back.
   But in some cases you want to change the default behaviour, e.g. declaring lower-case letters to be valid, as well.
   <br>
   **Optional task:** Implement the static member function `bool char_is_valid(char const)` in order to allow
   also 's' and 'w' as valid characters.
\endassignment
\solution
\snippet dna2_alphabet.cpp writable_alphabet
\endsolution

At this point the seqan3::alphabet concept should be modelled successfully and even seqan3::writable_alphabet
is fine because we implemented `assign_char`.
\snippet dna2_alphabet.cpp writable_alphabet_concept

## Shortcut: alphabet base template

Often it is not required to implement the entire class yourself, instead you can derive from seqan3::alphabet_base
which defines most things for you if you provide certain conversion tables. Read the documentation of
seqan3::alphabet_base for details. This implementation deduces `bool` as the smallest possible rank type.
\snippet dna2_derive_from_base.cpp dna2

# Further examples
## Implementation as enum class

This is an example of a minimal custom alphabet that provides implementations for all necessary customisation
points.

As an enum class the values already have an order and therefore the class models std::totally_ordered
without defining the (in)equality and comparison operators. Opposed to the examples above, we use only free functions
to model the functionality.
\snippet test/unit/alphabet/custom_alphabet_test.cpp my_alph

## Adaptation of a third party type {#howto_write_an_alphabet_custom}

This example is similar to the previous one, but assuming that you cannot add anything to the namespace of
the type that you wish to adapt.
In that case, you need to specialise the seqan3::custom::alphabet class template and provide the required functionality
as static members.

\snippet test/unit/alphabet/custom_alphabet3_test.cpp third_party_type

## Implementation of a non-default-constructible class

This is an example of a custom alphabet that is not default-constructible and that has a non-default overload for
seqan3::char_is_valid_for.

Please note that for the overloads of seqan3::alphabet_size and seqan3::char_is_valid_for our alphabet type has to
be wrapped into `std::type_identity<>` to be recognised by the customisation point objects, because our type does
not model std::is_nothrow_default_constructible after we have deleted the default constructor.

With the overload of seqan3::char_is_valid_for we allow assignment to the underlying boolean type
from '1', 't' and 'T' for value 1 as well as from '0', 'f' and 'F' for value 0.

\snippet test/unit/alphabet/custom_alphabet2_test.cpp my_alph

\note
You should really make your alphabet types [no-throw-default-constructible](\ref std::is_nothrow_default_constructible)
if you can!
