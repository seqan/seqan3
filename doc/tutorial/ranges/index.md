# Ranges {#tutorial_ranges}

[TOC]

This tutorial introduces the notion of *ranges*, a C++20 feature that SeqAn3 makes strong use of.

\tutorial_head{Moderate, 90 min, \ref tutorial_concepts,}

# Motivation

Traditionally most generic algorithms in the C++ standard library, like std::sort, take a pair of iterators
(e.g. the object returned by `begin()`).
If you want to sort a std::vector `v`, you have to call `std::sort(v.begin(), v.end())` and not `std::sort(v)`.
Why was this design with iterators chosen?
It is more flexible, because it allows e.g.:
  * sorting only all elements after the fifth one: `std::sort(v.begin() + 5, v.end())`
  * using non-standard iterators like reverse iterators: `std::sort(v.rbegin() + 5, v.rend())` (sorts in reverse order)

But this interface is less intuitive than just calling std::sort on the entity that you wish to sort and
it allows for more mistakes, e.g. mixing two incompatible iterators.
C++20 introduces the notion of *ranges* and provides algorithms that accept such in the namespace `std::ranges::`, e.g.
`std::ranges::sort(v)` now works if `v` is range – and vectors are ranges!

What about the two examples that suggest superiority of the iterator-based approach? In C++20 you can do the following:
  * sorting only all elements after the fifth one: `std::ranges::sort(std::views::drop(v, 5))`
  * sorting in reverse order: `std::ranges::sort(std::views::reverse(v))`

We will discuss later what `std::views::reverse(v)` does, for now it is enough to understand that it returns something
that appears like a container and that std::ranges::sort can sort it.
Later we will see that this approach offers even more flexibility than working with iterators.

# Ranges

*Ranges* are an abstraction of "a collection of items", or "something iterable". The most basic definition
requires only the existence of `begin()` and `end()` on the range.

There are different ways to classify ranges, one way is through the capabilities of its iterator.

## Range concepts

Ranges are typically either \link std::ranges::input_range input ranges \endlink (they can be read from) or
\link std::ranges::output_range output ranges \endlink (they can be written to) or both.
E.g. a `std::vector<int>` is both, but a `std::vector<int> const` would only be an input range.

\link std::ranges::input_range Input ranges \endlink have different *strengths* that are realised through more
refined concepts (i.e. types that model a stronger concept, always also model the weaker one):

| Concept                           | Description                                                 |
|-----------------------------------|-------------------------------------------------------------|
| std::ranges::input_range          | can be iterated from beginning to end **at least once**     |
| std::ranges::forward_range        | can be iterated from beginning to end **multiple times**    |
| std::ranges::bidirectional_range  | iterator can also move backwards with `--`                  |
| std::ranges::random_access_range  | you can jump to elements **in constant-time** `[]`          |
| std::ranges::contiguous_range     | elements are always stored consecutively in memory          |

For the well-known containers from the standard library this matrix shows which concepts they model:

|                                    | std::forward_list | std::list | std::deque | std::array | std::vector |
|------------------------------------|:-----------------:|:---------:|:----------:|:----------:|:-----------:|
| std::ranges::input_range           | ✅                | ✅        | ✅         | ✅         | ✅          |
| std::ranges::forward_range         | ✅                | ✅        | ✅         | ✅         | ✅          |
| std::ranges::bidirectional_range   |                   | ✅        | ✅         | ✅         | ✅          |
| std::ranges::random_access_range   |                   |           | ✅         | ✅         | ✅          |
| std::ranges::contiguous_range      |                   |           |            | ✅         | ✅          |

There are also range concepts that are independent of input or output or one of the above concept, e.g.
std::ranges::sized_range which requires that the size of a range can be computed and in constant time.

## Storage behaviour

**Containers** are the ranges most well known, they own their elements. SeqAn3 makes use of standard STL containers
like `std::vector`, but also implements some custom containers.

**Decorators** are ranges that are always defined on another range and decorate/annotate the underlying range
with additional information. They do not own the underlying range, but can contain member data of their own.

**Views** are ranges that are usually defined on another range and transform the underlying range
via some algorithm or operation.
Views do not own any data beyond their algorithm and the time it takes to construct, destruct or copy them should not
depend on the number of elements they represent. The algorithm is required to be lazy-evaluated so it is feasible to
combine multiple views. More on this below.

If you are confused about *decorators* vs *views*, think of decorators as "underlying range + data" and
views as "underlying range + algorithm".

The storage behaviour is orthogonal to the range concepts defined by the iterators mentioned above, i.e. you
can have a container that satisfies std::ranges::random_access_range (e.g. `std::vector` does, but `std::list`
does not) and you can have views or decorators that do so or don't. For some combinations of iterator capabilities
and storage behaviour there are extra concept definitions, e.g. seqan3::random_access_container.

# Views

As mentioned above, views are a specific kind of range.
They are incredibly useful and you will find them throughout the library.

## Lazy-evaluation

A key feature of views is that whatever transformation they apply, they do so at the moment you request an
element, not when the view is created.

\snippet doc/tutorial/ranges/range_snippets.cpp def

Here `v` is a view; creating it neither changes `vec`, nor does `v` store any elements.
The time it takes to construct `v` and its size in memory is independent of the size of `vec`.

\snippet doc/tutorial/ranges/range_snippets.cpp all

This will print "6", but the important thing is that resolving the first element of `v` to the last element of `vec`
happens **on-demand**.
This guarantees that views can be used as flexibly as iterators, but it also means that if the view performs an
expensive transformation, it will have to do so repeatedly if the same element is requested multiple times.


## Combinability

You may have wondered why we wrote

\snippet doc/tutorial/ranges/range_snippets.cpp rev_def

and not
```cpp
std::views::reverse v{vec};
```

That's because `std::views::reverse` is not the view itself, it's an *adaptor* that takes the underlying range
(in our case the vector) and returns a view object over the vector.
The exact type of this view is hidden behind the `auto` statement.
This has the advantage, that we don't need to worry about the template arguments of the view type, but more importantly
the adaptor has an additional feature: it can be *chained* with other adaptors!

\snippet doc/tutorial/ranges/range_snippets.cpp piped

What will this print?
\hint
It will print "4".
\endhint

In the above example the vector is "piped" (similar to the unix command line) into the reverse adaptor and then into
the drop adaptor and a combined view object is returned.
Note that accessing the 0th element of the view is still lazy, determining which element it maps to happens at the time
of access.

\assignment{Assignment 1: Fun with views I}
Look up the documentation of std::views::transform and std::views::filter.
Both take a invocable object as parameter, e.g. a lambda function.
std::views::transform applies the lambda on each element in the underlying range and std::views::filter
filter "removes" those elements that its lambda function evaluates to false for.

What does this imply for argument types and return types of the lambda functions?

\hint
The transform's lambda should return something of the same type as the input and the filter's lambda should return
true or false!
\endhint

Task: Create a view on `std::vector vec{1, 2, 3, 4, 5, 6};` that filters out all uneven numbers and squares the
remaining (even) values, i.e.
```cpp
std::vector vec{1, 2, 3, 4, 5, 6};
auto v = vec | // ...?

std::cout << *v.begin() << '\n'; // should print 4
```
\endassignment
\solution
\include doc/tutorial/ranges/range_solution1.cpp
\endsolution

## View concepts

Views are a specific kind of range that is formalised in the std::ranges::view concept.
Every view returned by a view adaptor models this concept, but which other range concepts are modeled by a view?

It depends on the underlying range and also the view itself.
With few exceptions, views don't model more/stronger range concepts than their underlying range (other than
std::ranges::view) and they try to preserve as much of the underlying range's concepts as possible.
For instance the view returned by `std::views::reverse` models std::ranges::random_access_range (and weaker concepts)
iff the underlying range also models the respective concept.
It never models std::ranges::contiguous_range, because the third element of the view is not located immediately after
the second in memory (but instead before the second).

Perhaps surprising to some, many views also model std::ranges::output_range if the underlying range does, i.e. **views
are not read-only**:

\snippet doc/tutorial/ranges/range_snippets.cpp assign_through

\assignment{Assignment 2: Fun with views II}
Have a look at the solution to the previous assignment (filter+transform).
Which of the following concepts do you think `v` models?

| Concept                          | yes/no? |
|----------------------------------|:-------:|
| std::ranges::input_range         |         |
| std::ranges::forward_range       |         |
| std::ranges::bidirectional_range |         |
| std::ranges::random_access_range |         |
| std::ranges::contiguous_range    |         |
|                                  |         |
| std::ranges::view                |         |
| std::ranges::sized_range         |         |
| std::ranges::output_range        |         |

\endassignment
\solution

| Concept                          | yes/no? |
|----------------------------------|:-------:|
| std::ranges::input_range         |   ✅    |
| std::ranges::forward_range       |   ✅    |
| std::ranges::bidirectional_range |   ✅    |
| std::ranges::random_access_range |         |
| std::ranges::contiguous_range    |         |
|                                  |         |
| std::ranges::view                |   ✅    |
| std::ranges::sized_range         |         |
| std::ranges::output_range        |         |

Surprised? Let's have a closer look at the std::views::filter view. The filter view only returns the value of the 
underlying range for which the given predicate evaluates to `true`. To know which value is an element of the filter
view, the view has to look at each of them. Thus, it must scan the underlying range value-by-value and cannot jump to an 
arbitrary location in constant time since it cannot know how many elements it had to skip without looking at them. 
Accordingly, the std::views::filter preserves only std::ranges::bidirectional_range, because it can scan the text in 
reverse order as well. Since the view cannot guarantee that the values lie in contiguous memory, it can also not 
preserve std::ranges::contiguous_range. Similarly, the view cannot model std::ranges::sized_range as it cannot determine
the number of values not filtered out in constant time.

The transform on the other hand produces a new element on every access (the result of the multiplication), therefore
`v` is not a std::ranges::output_range, you cannot assign values to its elements.
Note that this prevents modelling the std::ranges::contiguous_range as well because values are created on-demand and
are not stored in memory at all.
\endsolution

We provide overview tables for all our view adaptors that document which concepts are modelled by the views they return.

## Views in the standard library and in SeqAn

The standard library in C++20 provides a number of useful views and SeqAn provides many views, as well.
Most views provided by SeqAn3 are specific to biological operations, like seqan3::views::trim which trims sequences
based on the quality or seqan3::views::complement which generates the complement of a nucleotide sequence.
But SeqAn3 also provides some general purpose views.

Have a look at the \link views views-submodule \endlink to get an overview of SeqAn's views and also read through the
detailed description on that page now that you had a more gentle introduction.

\assignment{Assignment 3: Fun with views III}
Create a small program that
  1. reads a string from the command line (first argument to the program)
  2. "converts" the string to a range of seqan3::dna5 (Bonus: throw an exception if loss of information occurs)
  3. prints the string and it's reverse complement
  4. prints the six-frame translation of the string

Use views to implement steps 2.-4.
\endassignment
\solution
\include doc/tutorial/ranges/range_solution3.cpp
\endsolution

# Containers

containers are ranges that own their data.
SeqAn3 uses the standard library containers, like std::vector and std::list to store elements.
For certain use-cases we have introduced our own containers, though.

All standard library containers model std::ranges::forward_range (see above), but we have introduced container
concepts that encompass more of a containers interface.
Have a look at the API documentation of seqan3::container and unfold the inheritance diagram.
What can you learn about the different refinements and their relation to the range concepts?

## The bitcompressed vector

If you followed the alphabet tutorial closely, you will know that seqan3::dna4 needs only two bits to represent its state.
However, single objects are always at least a byte (eight bits) big in C++.
To store sequences of small alphabets more space-efficiently, we have developed seqan3::bitcompressed_vector.

Open the API documentation of seqan3::bitcompressed_vector, display the inheritance diagram and read through the
interface overview and the detailed description.

\assignment{Assignment 4: The bitcompressed vector}
Create a small program that asks the user for a size and then creates a vector of seqan3::dna4 of that size.
Add an argument parser flag that allows the user to decide whether std::vector or seqan3::bitcompressed_vector is used.
After creating the vector, print its size.

Measure and compare the amount of main memory that your program uses depending on the vector implementation.
On Linux based systems use `/usr/bin/time -v <program> <args>` and look for "Maximum resident set size".
(Not to be confused with the built-in Bash time command! So use the full path `/usr/bin/time`)

On macOS and BSD use `/usr/bin/time -l <program> <args>` and look for "maximum resident set size".

\note This command diplays the peak memory usage and only gives you a first impression. You can use [valgrind]
(http://valgrind.org/docs/manual/ms-manual.html) if you want to have a more detailed analysis of your memory
consumption.
\endassignment
\solution
\include doc/tutorial/ranges/range_solution4.cpp
\endsolution
