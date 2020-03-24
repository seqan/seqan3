# C++ Concepts {#tutorial_concepts}

[TOC]

This tutorial introduces "C++ Concepts", a feature of C++20 (and available to some extent in older GCC versions).
You will learn the terminology used in the context of concepts and how to use SeqAn's concepts in your application.

\tutorial_head{Moderate, 60 min, \ref setup\, \ref tutorial_argument_parser,}

This tutorial teaches the very basics of working with concepts. For more background and information on how to implement
your own concepts, we recommend:
  * A well-readable [paper](http://www.stroustrup.com/good_concepts.pdf) with motivation and historical background.
  * The (rather dense) [documentation on cppreference](https://en.cppreference.com/w/cpp/language/constraints).

# Constraints

## Motivation

One central design goal of SeqAn is to provide generic algorithms and data structures which can be used for different
types without reimplementing the same algorithms over and over again for particular types.
This has multiple benefits: improved maintainability due to an additional level of abstraction
and more importantly the ability to reuse the code with user provided types.
A familiar example for generic code is std::vector and the algorithms in the standard library.
They are *templates* which means that they can be *instantiated* with other types.
Most often the type cannot be arbitrary, because the template expects a particular interface from the type.

A SeqAn example is the local alignment algorithm.
It computes the best local match between two sequences over a finite alphabet.
The algorithm is generic in so far that it allows any alphabet that offers the minimal interface which
is used inside the algorithm (e.g. objects of the alphabet type must be equality comparable).
Before C++20, this could not be checked easily and using the interface with non-conforming types would result in
very hard to read compiler errors and consequently frustration of the user.
In the following part of the tutorial you will learn how to *constrain* such template arguments of generic functions
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
But of course it won't actually compile for all types, because the function template **implicitly requires** that
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

\attention The terse syntax in this form is not yet available in GCC7, GCC8 and GCC9.

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

\include doc/tutorial/concepts/overloading1.cpp

But as soon as we introduce another overload, the compiler will pick the "best" match:

\include doc/tutorial/concepts/overloading2.cpp

\assignment{Assignment 1: Static polymorphism with alphabets I}
Write a small program, similar to the one above with the following "skeleton":
```cpp
// which includes?

using seqan3::operator""_dna5;
using seqan3::operator""_aa27;

// Add one or more `void print` function template(s) here //

int main()
{
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
\include doc/tutorial/concepts/overloading_solution1.cpp
\endsolution

\assignment{Assignment 2: Static polymorphism with alphabets II}
Adapt your previous solution to handle nucleotides differently from the rest. For nucleotides, it should print both the value and its complement.
\endassignment
\solution
\include doc/tutorial/concepts/overloading_solution2.cpp
\endsolution

## Partial template specialisation

Similar to function template overloading it is possible to use concepts for partially specialising class and variable
templates.

\include doc/tutorial/concepts/specialisation.cpp

This is a typical example of a "type transformation trait".
It maps one type to another type; in this case it returns a type that is able to represent the square root of the
"input type".
This can be used in generic algorithms to hold data in different types depending on the type of the input –
in this case we could avoid half of the space consumption for unsigned integral types VS signed integral types.

\note The std::same_as used above is a concept with two template parameters.
It requires that both parameters are the same. The `static_assert` checks conditions at compile-time; it can be
used to verify whether a type or a combination of types model a concept. In the above case we can use the combination
to check the "return type" of the transformation trait.

# Concepts in SeqAn and this documentation

SeqAn uses concepts extensively, for specialisation/overloading, but also to prevent misuse of templates and to clearly
specify all public interfaces.
We prefer the intermediate syntax and additionally use the verbose expressions if necessary.
Unfortunately, doxygen, the system used to generate this documentation, does not handle C++ concepts very well, yet.
In some parts of the documentation concepts are called "interfaces", please don't let this confuse you.
And the "verbose syntax" introduced above is not visible at all in the automatically generated documentation.
That's why it's important to read the detailed documentation section where all requirements are documented.

Have a look at the documentation of seqan3::argument_parser::add_positional_option().
It has two template parameters, one seems unconstrained (`typename` in the signature) and one is constrained
(`validator` in the signature).
But in fact both are constrained as the detailed documentation reveals.

Now, follow the link to seqan3::validator. We will check in the next section whether you understand the
documentation for the concept.

# How to make your own type model a concept

## seqan3::validator

Remember the tutorial on \ref tutorial_argument_parser ? Let's implement our own validator that checks
if a numeric argument is an integral square (i.e. the user shall only be allowed to enter 0, 1, 4, 9...).

### Understanding the requirements

In the previous section you analysed seqan3::validator.
Do you understand the requirements formulated on that page?

\hint
In order to model the seqan3::validator, your custom validator must provide the following:

  1. It needs to expose a `value_type` type member which identifies the type of variable the validator works on.
     Currently, the SeqAn validators either have value_type `double` or `std::string`.
     Since the validator works on every type that has a common reference type to `value_type`, it enables a validator
     with `value_type = double` to work on all arithmetic values.
     \attention In order to be chainable, the validators need to share the same value_type!
  2. It has to be a [functor](https://stackoverflow.com/questions/356950/what-are-c-functors-and-their-uses), which
     basically means it must provide `operator()`.
  3. It has to have a member function `std::string get_help_page_message() const` that returns a string that can be
     displayed on the help page.

\endhint

### Formally satisfying the requirements

As we have noted previously, you can check if your type models seqan3::validator in the following way:

```cpp
struct custom_validator
{
    // ...
};

static_assert(seqan3::validator<custom_validator>);
```

To formally satisfy the requirements, your functions don't need the correct behaviour, yet.
Only the signatures need to be fully specified.

\assignment{Assignment 3: Custom validator I}
Implement enough of the above mentioned `struct custom_validator` for it to model seqan3::validator and pass
the check. You can use an empty `main()`-function for now.
\endassignment
\solution
\include doc/tutorial/concepts/custom_validator_solution1.cpp
\endsolution

### Implementing the functionality

The above implementation is of course not yet useful.
It should be usable with this main function:

\snippet doc/tutorial/concepts/custom_validator_solution2.cpp main

Try to think of the correct behaviour of this program.

It should print "Yeah!" for the arguments `-i 0`, `-i 4`, or `-i 144`; and/or `-j 0` or `-j 4`.

It should fail for the arguments `-i 3`; and/or `-j 144` or `-j 3`.

\assignment{Assignment 4: Custom validator II}
Implement your validator fully, i.e. make it throw seqan3::validation_error if the number provided is not a
square.
Also give a nice description for the help page.

\endassignment
\solution
\snippet doc/tutorial/concepts/custom_validator_solution2.cpp validator
\endsolution

You have now written your own type that is compatible with our constrained interfaces!
