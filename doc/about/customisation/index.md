# Customisation {#about_customisation}

<!-- SPDX-FileCopyrightText: 2006-2025 Knut Reinert & Freie Universität Berlin
     SPDX-FileCopyrightText: 2016-2025 Knut Reinert & MPI für molekulare Genetik
     SPDX-License-Identifier: CC-BY-4.0
-->

SeqAn clearly documents which of its interfaces are designed to be used with user-provided types and how
this can be realised.

[TOC]

# Generic interfaces

Templates are a core part of generic programming in C++, because they allow instantiating the code with different
types or constants. SeqAn's API consists mostly of templates for exactly this reason, although we try to make this
less visible by having the template parameters deduced in many scenarios.

We strive to document the exact requirements for every template parameter in the API documentation; often these
are also enforced via C++ concepts. Wherever other types meet the given requirements, you can use these
types as arguments – independent of whether they are standard library types, your own types or types from another
library. See \ref tutorial_concepts for more details on this.

Some requirements are more elaborate and depend on a specialisation or "overload" of another entity from the
SeqAn library. We call these "customisation points" (also known as "extension points").

# Customisation points

One such customisation point is seqan3::to_rank. Instead of the seqan3::semialphabet concept looking directly for
member functions or free functions in your namespace, it always checks for a valid definition of seqan3::to_rank which
in turn resolves to different possible implementations.

To customise one of our customisation points, follow the instructions in its API documentation. Typically, you have the
choice between the following options:

  1. Provide functionality as members.
  2. Provide friends or free functions inside the namespace of your type (will be picked up via
     [argument dependent lookup](https://en.cppreference.com/w/cpp/language/adl)).
  3. Specialise the respective template in `seqan3::custom::` for that type and provide the needed functionality as
     static members of that specialisation (this is an "upload space" for specialisations).

The priority is bottom to top (i.e. solutions to 3. will be preferred over 1. and 2.), but we strongly recommend to
use 1. and 2., because they are much easier to implement. Only use 3. if you adapt a third party's type and you cannot
add to that type's definition and/or namespace.

\warning
**Never** add anything (types, functions, variables...) to namespace `seqan3` and never explicitly specialise one
of our templates (except those in seqan3::custom) or overload one of our functions.

The \link howto_write_an_alphabet HowTo on creating your own alphabet \endlink provides many examples of how to
satisfy the requirements of customisation point objects.

More technical background on this topic can be found here:

  1. [Original blog article](https://ericniebler.com/2014/10/21/customization-point-design-in-c11-and-beyond/) with
     background and the first "Niebloid" design.
  2. [A more recent article](https://quuxplusone.github.io/blog/2018/03/19/customization-points-for-functions/) that
     also introduces the notion of "upload" namespaces.

There is also a proposed [language extension](http://www.open-std.org/jtc1/sc22/wg21/docs/papers/2018/p1292r0.html)
to handle customisation points in a future version of C++.
