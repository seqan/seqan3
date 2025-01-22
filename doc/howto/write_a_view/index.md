# How to write a view {#howto_write_a_view}

<!-- SPDX-FileCopyrightText: 2006-2025 Knut Reinert & Freie Universität Berlin
     SPDX-FileCopyrightText: 2016-2025 Knut Reinert & MPI für molekulare Genetik
     SPDX-License-Identifier: CC-BY-4.0
-->

[TOC]

This HowTo documents how to write a view using the standard library and some helpers from SeqAn.

\tutorial_head{Difficult, 120 min, \ref tutorial_concepts\,\ref tutorial_ranges, }

\note
Some of the links from this HowTo only resolve in the developer documentation because
they refer to entities from the seqan3::detail namespace.
We recommend you open this tutorial from [the developer documentation](https://docs.seqan.de/seqan3/main_dev/).

# Motivation

We have introduced "views" in \ref tutorial_ranges.
You can do many things with the views provided by the standard library and those shipped with SeqAn, but in
certain situations you will want to define your own view.
This page will teach you the basics of defining your own view.

# What makes a view?

A view is a type of `std::ranges::range` that also models `std::ranges::view`.
The additional requirements of `std::ranges::view` can be vaguely summarised as "not holding any own data" or at least
not holding data that is relative in size to the number of elements in the view (e.g. a vector cannot be a view,
because its size in memory depends on the number of elements it represents).

A simple example of a view is `std::ranges::subrange`. It can be constructed from a pair of iterators or more precisely
an *iterator* and a *sentinel*. `begin()` always returns an iterator, but the type returned by `end()` (the "sentinel")
does not need to be of the same type as the iterator – as long as they are comparable.
This view then holds exactly the iterator-sentinel-pair as its state and nothing else.

But `std::ranges::subrange` does not yet facilitate any "composing-behaviour" that you have seen in
\ref tutorial_ranges tutorial,
`std::ranges::subrange` is simply a type that can be constructed from an iterator-sentinel-pair, you cannot "pipe"
anything into it.
Most views are adaptors on other views, e.g. `std::ranges::transform_view` wraps an existing view and applies
an element-wise transformation on-demand.
You can directly construct `std::ranges::transform_view` from another view or from non-view ranges that can be wrapped
into a view (e.g. references to containers):

```cpp
std::vector vec{1, 2, 3, 4, 5};

auto l = [] (int i) { return i + 1; };

std::ranges::transform_view v{vec, l};
```

But this syntax gets difficult to read when you create "a view(from a view(from a view()))".
That's why for every view that adapts an existing view we have an additional *adaptor object*, usually available
in a `views::` sub-namespace:

```cpp
std::vector vec{1, 2, 3, 4, 5};

auto l = [] (int i) { return i + 1; };

auto v = vec | std::views::transform(l);
```

This adaptor object (`std::views::transform`) provides the pipe operator and returns an object of the actual view type
(`std::ranges::transform_view`).
The pipe operator allows us to chain multiple adaptors similar to the unix command line.
We will discuss the details of these adaptor objects in the following sections.

# Custom range adaptor objects

Read [section 24.7 and 24.7.1 of the C++ standard](https://eel.is/c++draft/range.adaptors).

The wording of the standard needs some getting used to, but some important notes for us are:
  1. You can pipe a viewable range into series of *adaptor closure objects* and will get back a view (this is also
  what we did above):
```cpp
std::vector vec{1, 2, 3, 4, 5};
auto v = vec | std::views::transform([] (int i) { return i + 1; })
             | std::views::filter([] (int i) { return i % 2 == 0; });
// v is a view, you can iterate over it!
```
  2. The adaptor objects support function-style usage, too, although it only reduces readability:
```cpp
std::vector vec{1, 2, 3, 4, 5};
auto v = std::views::filter(std::views::transform(vec, [] (int i) { return i + 1; }), [] (int i) { return i % 2 == 0; });
// v is a view, you can iterate over it!
```

  3. You can create a new *adaptor closure object* from an adaptor object that requires parameters by providing those:
```cpp
std::vector vec{1, 2, 3, 4, 5};
auto a = std::views::transform([] (int i) { return i + 1; });
// a is an adaptor and can be used as such:
auto v = vec | a;
```

  4. You can create a new *adaptor closure object* from existing adaptor closure objects via `|` but without providing
  a range:
```cpp
std::vector vec{1, 2, 3, 4, 5};
auto a = std::views::transform([] (int i) { return i + 1; })
       | std::views::filter([] (int i) { return i % 2 == 0; });
// a is an adaptor and can be used as such:
auto v = vec | a;
```

| Terminology                    | Description                                                               | Example                                                                                                           |
|--------------------------------|---------------------------------------------------------------------------|-------------------------------------------------------------------------------------------------------------------|
| "view"                         | A type that models std::ranges::view.                                     | `std::ranges::filter_view<std::subrange<int const *, int const *>>`                                               |
| "view adaptor object"          | Creates a view; can be combined with other adaptors.                      | `std::views::reverse` and <br> `std::views::filter` and <br> `(std::views::filter([] (int) { return true; }))`       |
| "view adaptor closure object"  | A "view adaptor object" that requires no paramaters other than the range. | `std::views::reverse` and \strike{<tt>std::views::filter</tt> and} `(std::views::filter([] (int) { return true; }))` |

In many cases where you are planning on creating "a new view", it will be sufficient to use the previously
mentioned techniques to just create "a new adaptor object" and not having to specify the actual view type yourself.

Let's look at some examples!

\assignment{Exercise 1: Your first custom adaptor object}
In the alphabet module you learned that ranges of alphabets are not implicitly
convertible to `char` so you **cannot** print a `std::vector<seqan3::dna5>` via `std::cout`¹.
You also know that you can call `seqan3::to_char` on every object that models `seqan3::alphabet` which will convert it
to `char` or a similar type.

We want to do the following:

\snippet doc/howto/write_a_view/view_exercise1.cpp start
```cpp
// your implementation goes here
```
\snippet doc/howto/write_a_view/view_exercise1.cpp end

Define a *range adaptor object* using an existing adaptor which applies a concrete transformation
(calling seqan3::to_char) on every element.

\hint
You need to need use `std::views::transform` and you need to set a fixed transformation function. `std::views::transform`
takes an object that models `std::regular_invocable`, e.g. a lambda function with empty capture `[]`.
\endhint

<small>¹ You *can* print via `seqan3::debug_stream`, but let's ignore that for now. </small>
\endassignment

\solution

\include doc/howto/write_a_view/view_exercise1.cpp

You simply define your adaptor type as `auto` and make it behave like `std::views::transform`, except that you
"hard-code" the lambda function that is applied to each element.
Since your adaptor object now takes a range as the only parameter, it is an adaptor closure object.

The object is marked as `constexpr` because the adaptor object itself never changes, it only
provides `operator()` and `operator|` that each return a specialisation of `std::ranges::transformation_view`.
\endsolution

\assignment{Exercise 2: Combining two existing adaptor objects}

Study the `seqan3::nucleotide_alphabet`. It states that you can call `seqan3::complement` on all nucleotides
which will give you <tt>'A'_dna5</tt> for <tt>'T'_dna5</tt> a.s.o. Think about how you can adapt the previous solution
to write a view that transforms ranges of nucleotides into their complement.

BUT, we are also interested in *reversing* the range which is possible with `std::views::reverse`:

\snippet doc/howto/write_a_view/view_exercise2.cpp start
```cpp
// your implementation goes here
```
\snippet doc/howto/write_a_view/view_exercise2.cpp end

Define a *range adaptor object* that presents a view of the reverse complement of whatever you pipe into it.
\endassignment

\solution

\include doc/howto/write_a_view/view_exercise2.cpp

The adaptor consists of `std::views::reverse` combined with `std::views::transform`. This time the lambda just
performs the call to `seqan3::complement`.
\endsolution

# A full custom view implementation

Using existing adaptors only works to a certain degree, sometimes you will need to implement a full view.
As we have seen above, it is easy to implement a view that presents the complement of a nucleotide range just using
existing adaptors.
For simplicity, we will still use this as an example and implement a non-generic transform view in
the next steps.

A full view implementation typically consists of the following components:

  1. The view class template that takes the underlying range as a template parameter.
     * The class should be derived from `std::ranges::view_interface` which will take care of a lot of boilerplate for
       us.
     * The actual functionality of the view is implemented in its iterator and sentinel.
  2. The *adaptor object* which allows for the piping behaviour.

## 1. The view class template

A good thing to start with is to think about the iterator and the sentinel of your view.
The iterator and/or its relation to the sentinel is the way we implement the behaviour that is specific to
our view.
We will start with implementing the iterator separately and later integrate it into the view.

### Iterator

Since we know that in our current usecase we have exactly one element in our view for every element in the
underlying range, we don't need to change the relation between iterator and sentinel, i.e. "begin == end" on
our view iff "begin == end" on the underlying range.
This indicates that we can use the sentinel of the underlying range as-is.

\assignment{Exercise 3: Your first iterator}

Have a peak at the \ref tutorial_concepts tutorial again and study the `std::forward_iterator` concept thoroughly.
You will now have to implement your own forward iterator.

\snippet doc/howto/write_a_view/solution_iterator.cpp start

```cpp
    // YOUR IMPLEMENTATION GOES HERE (OPERATORS, MEMBER TYPES...)
```
\snippet doc/howto/write_a_view/solution_iterator.cpp end

In order to re-use functionality of the underlying range's iterator type you can inherit from it
(`std::ranges::iterator_t` returns the iterator type).¹ In the end, you should be able
to iterate over the underlying range and print elements like you would with the original iterator, i.e. your iterator
shall behave exactly as the original in this regard (no transformation, yet).

Some things to keep in mind for the implementation:
  * When defining a type template (like above), member **types** are not inherited implicitly (member functions are).
  * Another reason "just inheriting" is not sufficient, is that some inherited functions return objects or references
    to the base type, not your type (which is required by the `std::forward_iterator` concept).
  * If you choose to wrap the underlying iterator instead of inheriting, you will need to define and "forward"
    a few more member functions, but it's good practice and will help you understand the concept.
  * The `static_assert` will just return `true` or `false`, to get more detailed information on why your type does
    not (yet) model the concept, implement a constrained function template and pass an object of your type to that
    template.
  * Your iterator needs to be comparable to the sentinel of the underlying range, you can access that type via
    `std::ranges::sentinel_t<urng_t>`.


<small>¹ In some situations it might be better to wrap the underlying iterator instead of inheriting it, i.e.
save a copy of the underlying iterator as a data member of your iterator.
A reason could be that you *don't* want to inherit some members or want to prevent implicit convertibility.</small>

\endassignment

\solution
\snippet doc/howto/write_a_view/solution_iterator.cpp start
\snippet doc/howto/write_a_view/solution_iterator.cpp solution1a
```cpp
    using reference             = typename std::iterator_traits<base_t>::reference;

```
\snippet doc/howto/write_a_view/solution_iterator.cpp solution1b
\snippet doc/howto/write_a_view/solution_iterator.cpp end

The program prints "G A T T A C A ".
\endsolution

\assignment{Exercise 4: A transforming iterator}

In the previous assigment you have created a working – but pointless – iterator.
It does not do anything differently from the original.

Your task now is to implement the "complementing" behaviour, i.e.
  * your iterator should only work on ranges of `seqan3::nucleotide_alphabet`
  * when accessing an element through the iterator it should not return the element from the underlying range but
    instead its complement

Think about which operator is responsible for returning the element and be careful with the return type
of that operator as your transformation might make a change necessary.
\endassignment

\solution

1. The restriction on the alphabet type is done via a `static_assert`:
\snippet doc/howto/write_a_view/solution_iterator.cpp static_assert
You could have done this via an additional template constraint, too, but `static_assert` gives you the opportunity
to give a readable message in case of an error.¹

2. The operator that needs to call the `seqan3::complement` function is `operator*`:
\snippet doc/howto/write_a_view/solution_iterator.cpp dereference

3. As previously noted, care needs to be taken with this function's return type:
\snippet doc/howto/write_a_view/solution_iterator.cpp reference


Here is the full solution:
\hint
\include doc/howto/write_a_view/solution_iterator.cpp

The program prints "C T A A T G T "

\endhint

<small>¹ This is only recommended when you do *not* want to allow a different specialisation of the template to
cover the excluded case.</small>

\endsolution

You now have a working iterator, although it still lacks the capabilities of `std::bidirectional_iterator`,
`std::random_access_iterator` and `std::contiguous_iterator`.
When designing views, you should always strive to preserve as much of the capabilities of the underlying range
as possible.

Which of the mentioned concepts do you think your iterator could be made to implement? Have a look at the respective
documentation.

\hint
It could be designed to be a `std::random_access_iterator` (and thus also `std::bidirectional_iterator`) when the
underlying range is, because jumping on your iterator/view can be done in
constant time iff it can be done in constant time on the underlying range (you just jump to the n-th element
of the underlying range and perform your transformation one that).

However, it can never model `std::contiguous_iterator` because that would imply that the elements are adjacent to each
other in memory (the elements of our view are created on demand and are not stored in memory).
\endhint

If you have looked at the `std::random_access_iterator`, you will have seen that it is quite a bit of work to implement
all the operators, many of whom just need to be overloaded to fix the return type.
To make this a little bit easier SeqAn provides `seqan3::detail::inherited_iterator_base`, it fixes the issue with the
return type via CRTP.
A solution to the previous exercise looks like this:

\snippet doc/howto/write_a_view/solution_view.cpp iterator
\snippet doc/howto/write_a_view/solution_view.cpp main_it
\snippet doc/howto/write_a_view/solution_view.cpp end

### The view class

We now implement the view in several steps:

\snippet doc/howto/write_a_view/solution_view.cpp view_header

Like the iterator, the view is derived from a CRTP base class that takes care of defining many members
for us, e.g. `.size()`, `.operator[]` and a few others.

\snippet doc/howto/write_a_view/solution_view.cpp view_private

The only data member the class holds is a copy of the underlying range.
As you may have noted above, our class only takes underlying ranges that model
std::ranges::view.
This might seem strange; after all we want to apply the view to a vector of which we know that it
is not a view, but we will clear this up later.

\snippet doc/howto/write_a_view/solution_view.cpp view_member_types

The only member types that we define here are the definitions of the iterators which are just the iterator
we have defined before.
Note that we have `const_iterator` in addition to `iterator` which const-qualified member functions return,
because in a const-context the `urange` data member will be const so we cannot return the mutable `iterator` from it.

Many ranges like the standard library containers also present the member types of the iterator, i.e. `value_type`,
`reference` a.s.o, but this is not required to model any of the range concepts.

\snippet doc/howto/write_a_view/solution_view.cpp view_begin

These functions are the same member functions you know from `std::vector`, they return objects of the previously
defined iterator types that are initialised with the begin iterator from the underlying range.

\snippet doc/howto/write_a_view/solution_view.cpp view_end

The implementation for `end()` is similar except that for our range the sentinel type (the return type of `end()`) is
the same as of the underlying range, we just pass it through.

For many more complex views you will have to define the sentinel type yourself or derive it from the underlying
type in a similar manner to how we derived the iterator type.
Often you can use `std::default_sentinel_t` as the type for your sentinel and implement the
"end-condition" in the iterator's equality comparison operator against that type.

\snippet doc/howto/write_a_view/solution_view.cpp view_constructors

We have two constructors, one that takes an the underlying type by copy and moves it into the data member
(remember that since it is a view, it will not be expensive to copy – if it is copied).

The second constructor is more interesting, it takes a `std::ranges::viewable_range` which is defined as
being either a `std::ranges::view` or a reference to `std::ranges::range` that is not a view
(e.g. `std::vector<char> &`).
Since we have a constructor for `std::ranges::view` already, this one explicitly handles the second case
and goes through `std::views::all` which wraps the reference in a thin view-layer.
Storing only a view member guarantess that our type itself is also cheap to copy among other things.

Note that both of these constructors seem like generic functions, but they just handle the underlying type or a
type that turns into the underlying when wrapped in `std::views::all`.

\snippet doc/howto/write_a_view/solution_view.cpp view_deduction_guide

To easily use the second constructor we need to provide a type deduction guide.

Here is the full solution:
\hint
\snippet doc/howto/write_a_view/solution_view.cpp iterator
\snippet doc/howto/write_a_view/solution_view.cpp view_header
\snippet doc/howto/write_a_view/solution_view.cpp view_private
\snippet doc/howto/write_a_view/solution_view.cpp view_member_types
\snippet doc/howto/write_a_view/solution_view.cpp view_constructors
\snippet doc/howto/write_a_view/solution_view.cpp view_begin
\snippet doc/howto/write_a_view/solution_view.cpp view_end
```cpp
};
```
\snippet doc/howto/write_a_view/solution_view.cpp view_deduction_guide
\snippet doc/howto/write_a_view/solution_view.cpp main_it
\snippet doc/howto/write_a_view/solution_view.cpp main_range
\snippet doc/howto/write_a_view/solution_view.cpp end

The program prints
```
C T A A T G T
CTAATGT
```
\endhint

## 2. The adaptor object

The *adaptor object* is a *function object* also called functor. This means we define a type with the respective
operators and then create a global instance of that type which can be used to invoke the actual functionality.

### Adaptor type definition

The adaptor has the primary purpose of facilitating the piping behaviour, but it shall also allow for
function/constructor-style creation of view objects, therefore it defines two operators:

\snippet doc/howto/write_a_view/solution_view.cpp adaptor_type_definition

The first operator is very straight-forward, it simply delegates to the constructor of our view so that
`views::my(FOO)` is identical to `my_view{FOO}`.

The second operator is declared as friend, because the left-hand-side of the `operator|` is generic, it's how
the range is handled in snippets like `auto v = vec | views::my`.

\note
This adaptor type does not yet provide the ability to combine multiple adaptors into a new adaptor, it only
handles ranges as left-hand-side input to `operator|`.

Our example adaptor type definition is rather simple, but for views/adaptors that take more parameters it gets quite
complicated quickly. Therefore SeqAn provides some convenience templates for you:
```cpp
// in our example, this is all you need:
//                                      your view type goes here ↓
using my_view_fn = seqan3::detail::adaptor_for_view_without_args<my_view>;
```

See `seqan3::detail::adaptor_base`, `seqan3::detail::adaptor_for_view_without_args` and
`seqan3::detail::adaptor_from_functor` for more details.

### Adaptor object definition

The adaptor object is simply an instance of the previously defined type:

\snippet doc/howto/write_a_view/solution_view.cpp adaptor_object_definition

As noted above, we place this object in a `views::` sub-namespace by convention.
Since the object holds no state, we mark it as `constexpr` and since it's a global variable we also mark it as
`inline` to prevent linkage issues.

Finally we can use our view with pipes and combine it with others:

\snippet doc/howto/write_a_view/solution_view.cpp main_adaptor

Here is the full, final solution:
\hint
\include doc/howto/write_a_view/solution_view.cpp

The program prints:
```
C T A A T G T
CTAATGT
TGTAATC
```

\endhint
