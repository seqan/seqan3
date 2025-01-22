# C++ Concepts {#tutorial_concepts}

<!-- SPDX-FileCopyrightText: 2006-2025 Knut Reinert & Freie Universität Berlin
     SPDX-FileCopyrightText: 2016-2025 Knut Reinert & MPI für molekulare Genetik
     SPDX-License-Identifier: CC-BY-4.0
-->

***Learning Objective:***

This tutorial introduces "C++ Concepts", a feature of C++20 (and available to some extent in older GCC versions).
You will learn the terminology used in the context of concepts and how to use SeqAn's concepts in your application.

\tutorial_head{Moderate, 60 min,
               \ref setup\,
               \ref tutorial_argument_parser,
               [Concepts (cppreference)](https://en.cppreference.com/w/cpp/language/constraints)}

This tutorial teaches the very basics of working with concepts. For more background and information on how to implement
your own concepts, we recommend:
  * A well-readable [paper](https://www.stroustrup.com/good_concepts.pdf) with motivation and historical background.
  * The (rather dense) [documentation on cppreference](https://en.cppreference.com/w/cpp/language/constraints).

[TOC]

# Constraints

## Motivation

One central design goal of SeqAn is to provide generic algorithms and data structures that can be used for different
types without reimplementing the same algorithms over and over again for particular types.
This has multiple benefits: improved maintainability due to an additional level of abstraction
and, more importantly, the ability to reuse the code with user provided types.
A familiar example for generic code is std::vector and the algorithms in the standard library.
They are *templates* which means that they can be *instantiated* with other types.
Most often the type cannot be arbitrary, because the template expects a particular interface from the type.

A SeqAn example is the local alignment algorithm.
It computes the best local match between two sequences over a finite alphabet.
The algorithm is generic in so far that it allows any alphabet that offers the minimal interface which
is used inside the algorithm (e.g. objects of the alphabet type must be equality comparable).
Before C++20, this could not be checked easily and using the interface with non-conforming types would result in
very hard to read compiler errors and consequently frustration of the user.
In the following part of the tutorial, you will learn how to *constrain* such template arguments of generic functions
and data structures and how this can have a huge impact on your code.

Here's a shorter example:

```cpp
template <typename t>
t add(t const v1, t const v2)
{
    return v1 + v2;
}

int main()
{
    return add(1, 3); // instantiates add<int>()
}
```

The template parameter `t` is said to be *unconstrained*, in theory it can be instantiated with any type.
But of course, it won't actually compile for all types because the function template **implicitly requires** that
types provide a `+` operator.
If a type is used that does not have a `+` operator, this implicitness causes the compiler to fail at the place where
such operator is used – and not at the place the template is instantiated.
This leads to very complex error messages for deeply nested code.

*Constraints* are a way of making requirements of template arguments **explicit**.
Constraints can be formulated ad-hoc, but this tutorial only covers *concepts*.
The interested reader can check the [documentation](https://en.cppreference.com/w/cpp/language/constraints) to learn
about ad-hoc definitions.
Concepts are a set of constraints with a given name.
Let's assume there is a concept called `Addable` that requires the existence of a `+` operator (as previously mentioned
the syntax for defining concepts is not covered here).
The following snippet demonstrates how we can constrain our function template, i.e. make the template immediately
reject any types that don't satisfy the requirement:

```cpp
template <Addable t>
t add(t const v1, t const v2)
{
    return v1 + v2;
}

int main()
{
    return add(1, 3); // instantiates add<int>()
}
```

The only difference is that we have replaced `typename` with `Addable`.
If you plug in a type that does not model `Addable`, you will get a message stating exactly that and not a cryptic
template backtrace.

The standard library provides a set of [predefined concepts](https://en.cppreference.com/w/cpp/concepts).
For our example above, the std::integral concept could have been used.

## Syntax variants

Depending on the complexity of your constraint statements, three different syntaxes are available to enforce
constraints; all of the following are equivalent.

(1) The "verbose syntax", especially useful when enforcing multiple constraints:

```cpp
template <typename t1, typename t2>
    requires std::integral<t1> && std::integral<t2> // && MyOtherConcept<t1>
auto add(t1 const v1, t2 const v2)
{
    return v1 + v2;
}
```

(2) The "intermediate syntax":
```cpp
template <std::integral t1, std::integral t2>                       // one constraint per type
auto add(t1 const v1, t2 const v2)
{
    return v1 + v2;
}
```

(3) The "terse syntax":
```cpp
auto add(std::integral auto const v1, std::integral auto const v2)  // one constraint per type
{
    return v1 + v2;
}
```

Different constraints can be applied to different template parameters and a single template parameter can be constrained
by multiple concepts.
Syntaxes can also be combined:
```cpp
template <std::integral t1, std::integral t2>
    // requires MyOtherConcept<t1>
auto add(t1 const v1, t2 const v2)
{
    return v1 + v2;
}
```

# Terminology

  * Template arguments can be ***constrained***.
  * A named set of constraints is a ***concept***.
  * A type that satisfies all requirements of a concept is said to ***model*** said concept.
  * A *concept* that is composed of another concept and additional constraints is said to ***refine*** said concept(s).

Some people confuse concepts with *interfaces*.
Both can be used as an abstraction of concrete types, but interfaces have to be inherited from. → the abstraction
is explicit in the definition of the type.
Concepts on the other hand "describe properties from the outside". → types don't need to be related and don't need
to "know about the concept" to model it.

Furthermore, the polymorphism possible with concepts (see below) is faster, because it is resolved at compile-time while
interface inheritance is resolved at run-time.

# Overloading and specialisation

In generic programming, "function overloading" and "template specialisation" play an important role.
They allow providing generic interfaces and (gradually) more specialised implementations for specific types or groups
of types.

## Function (template) overloading

When a function is overloaded and multiple overloads are valid for a given/deduced template argument, the
*most-refined* overload is chosen:

\include doc/tutorial/03_concepts/overloading1.cpp

But as soon as we introduce another overload, the compiler will pick the "best" match:

\include doc/tutorial/03_concepts/overloading2.cpp

\assignment{Assignment 1: Static polymorphism with alphabets I}
Write a small program, similar to the one above with the following "skeleton":
```cpp
// which includes?

// Add one or more `void print` function template(s) here //

int main()
{
    using namespace seqan3::literals;

    auto d = 'A'_dna5;
    auto a = 'L'_aa27;
    auto g = seqan3::gap{};

    print(d);
    print(a);
    print(g);
}
```

The `print` function (template) should print for every object `v` passed to it the result of `to_char(v)` and it should
be constrained to only accepts types that model seqan3::alphabet.
Try calling `print` with a different type, e.g. `int` to make sure that it does.
\endassignment
\solution
\include doc/tutorial/03_concepts/overloading_solution1.cpp
\endsolution

\assignment{Assignment 2: Static polymorphism with alphabets II}
Adapt your previous solution to handle nucleotides differently from the rest. For nucleotides, it should print both the value and its complement.
\endassignment
\solution
\include doc/tutorial/03_concepts/overloading_solution2.cpp
\endsolution

## Partial template specialisation

Similar to function template overloading it is possible to use concepts for partially specialising class and variable
templates.

\include doc/tutorial/03_concepts/specialisation.cpp

This is a typical example of a "type transformation trait".
It maps one type to another type; in this case, it returns a type that is able to represent the square root of the
"input type".
This can be used in generic algorithms to hold data in different types depending on the type of the input –
in this case, we could avoid half of the space consumption for unsigned integral types VS signed integral types.

\note The std::same_as used above is a concept with two template parameters.
It requires that both parameters are the same. The `static_assert` checks conditions at compile-time; it can be
used to verify whether a type or a combination of types model a concept. In the above case, we can use the combination
to check the "return type" of the transformation trait.

# Concepts in SeqAn and this documentation

SeqAn uses concepts extensively, for template specialisation/overloading, to avoid misuse and improve error messages.
Unfortunately, doxygen, the system used to generate this documentation, does not handle C++ concepts very well, yet.
That's why it's important to read the detailed documentation section of the constrained type, where we try to document
the requirements manually.
In some parts of the documentation concepts are called "interfaces", please don't let this confuse you.

<!-- To prevent misuse of templates and to clearly specify all public interfaces we use the concepts within
`static_assert`s in order to provoke more readable error messages. The thereby enforced requirements are also manually
documented with the respective instances. WE DO NOT ACTUALLY DO THIS... -->

## Example: seqan3::bitpacked_sequence

The class `seqan3::bitpacked_sequence<alphabet_type>` behaves just like `std::vector<alphabet_type>` but has an internal representation where multiple
values are packed into a single byte/word to save space. Also analog to `std::vector`, not every `alphabet_type` can
be used. To avoid misuse and weird error messages, the type is constrained.

Have a look at the documentation of [`seqan3::bitpacked_sequence`](http://docs.seqan.de/seqan3/main_user/classseqan3_1_1bitpacked__sequence.html).
It has one constrained template parameter.
Do you understand the requirements imposed on `alphabet_type` when using the
[`seqan3::bitpacked_sequence`](http://docs.seqan.de/seqan3/main_user/classseqan3_1_1bitpacked__sequence.html)?

\hint
In order to use the `seqan3::bitpacked_sequence` the  `alphabet_type` must model the following:

  1. It needs to model [`std::regular`](https://en.cppreference.com/w/cpp/concepts/regular), a stl concept.
     This only enforcing two other concepts: `std::semiregular<T> && std::equality_comparable<T>`.
     * `std::semiregular<T>` makes sure that your type is default initialisable (e.g. `int i{};`).
     * `std::equality_comparable<T>` makes sure you can compare your type with `==` (e.g. `i == j`).

     It makes sense that in order to save a range of letters (of type `alphabet_type`), you need them to be
     default initialisable, for example s.t. you can easily resize your container.
     Additionally, `seqan3::bitpacked_sequence` needs the `alphabet_type` to be comparable, in order be equality
     comparable itself (e.g. you can do `bit_seq_1 == bit_seq_2`).

  2. It needs to model [`seqan3::writable_semialphabet`], a seqan3 concept.
     This again enforces two things:
     * `seqan3::assign_rank_to` needs to be defined for objects of this type.
     * the type shall model `seqan3::semialphabet`,
       which in summary enforces that your type is ordered (comparable via `<`), shall be efficiently copyable and
       you should be able to call `seqan3::alphabet_size(c)` and `seqan3::to_rank(c)` (assuming `c` is of type `alphabet_type`).

\endhint

Of course, all seqan3 alphabets model the requirements and can be used with the `seqan3::bitpacked_sequence`.

But what happens if a type you would like to use does not model `seqan3::writable_semialphabet` (because obviously
this concept is very SeqAn specific)?

You can learn how to make your own alphabet model the SeqAn requirements in \ref howto_write_an_alphabet

In order to understand what "make a type model a concept" means in practical terms, let's look at an easier
example in the next section.

# Satisfying a concept

Let's say you have the following concept called `fooger`:

\snippet doc/tutorial/03_concepts/model_a_concept.cpp concept

Do you understand the requirements?

\hint
  1. The type `T` needs to model `has_foo<T>`
     Which again has two requirements:
     requirement 1: The type `T` has to have a *type member* called `FOO`
     requirement 2: The type `T` has to have a *member variable* calles `foo`
  2. `std::same_as` is a concept that checks whether two types are exaclty the same.
     Thus, `fooger` requires, that the *type member* `T::FOO` is `int`.
\endhint

\assignment{Assignment 4: Make a type model a concept}

Copy over the concept into a new `.cpp` file.

Add a type `my_type` that models the requirements, s.t.

\snippet doc/tutorial/03_concepts/model_a_concept.cpp main

prints `1`.

Hint: Don't forget to include the `seqan3::debug_stream` via `#include <seqan3/core/debug_stream.hpp>`.

\endassignment
\solution
\include doc/tutorial/03_concepts/model_a_concept.cpp
\endsolution
