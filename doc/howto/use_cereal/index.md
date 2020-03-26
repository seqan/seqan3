# How to serialise SeqAn data structures {#howto_use_cereal}

[TOC]

This HowTo shows how to serialise data structures with [cereal](https://uscilab.github.io/cereal/). Every SeqAn data
structure which is is marked as cerealisable can be used.

\tutorial_head{Easy, 15 min, No prerequisites , }

# Motivation

Storing and loading data, for example a [FM-Index](https://en.wikipedia.org/wiki/FM-index), is a common use case. Thanks
 to the cereal library doing so is incredible easy.  This page will show you how to use cereal in SeqAn. As an example
data structure we will use a `std::vector`, but as already mentioned, any SeqAn data structure that is documented
as cerealisable can be used.

# Storing

Storing a data structure is as easy as using the `cereal::BinaryOutputArchive`.
In order to use it, you need to include

\snippet doc/howto/use_cereal/store.hpp binary_include

Note that stl data types, like `std::vector`, are not automatically usable with cereal. You need to include the
respective header, e.g.
\snippet doc/howto/use_cereal/store.hpp vector_include

In SeqAn, our data structures have integrated cereal support so there is no need to include an extra header.

\include doc/howto/use_cereal/store.hpp

# Loading

Loading a data structure is as easy as using the `cereal::BinaryInputArchive`.
In order to use it, you need to include

\snippet doc/howto/use_cereal/load.hpp binary_include


\include doc/howto/use_cereal/load.hpp

# Storing & Loading in the same function

In the example above loading and storing was encapsulated in separated functions. It is possible to use
`cereal::BinaryInputArchive` and `cereal::BinaryOutputArchive` in one function, but then it is necessary to encapsulate
each in an individual scope (using extra braces`{}`). The reason for this is that the output/input stream handle of an
archive is closed on deconstruction and only then is the file properly written an accessible by another filehandle.
