# Customisation {#about_customisation}

SeqAn clearly documents which of its interfaces are designed to be used with user-provided types and how
this can be realised.

[TOC]

# Generic interfaces

Templates are a core part of generic programming in C++, because they allow instantiating the code with different
types or constants. SeqAn's API consists mostly of templates for exactly this reason, although we try to make this
less visible by having the template parameters deduced in many scenarios.

We strive to document the exact requirements for every template parameter in the API documentation; often these
are also enforced via C++ concepts. Wherever other types meet the given requirements, you can use these
type as arguments – independent of whether they are standard library types, your own types or types from another
library. See \ref tutorial_concepts for more details on this.

Some requirements are more elaborate and depend on a specialisation or "overload" of another entity from the
SeqAn library. We call these "customisation points" (also known as "extension points").

# Customisation points

One such customisation point is seqan3::to_rank. Instead of the seqan3::Semialphabet concept looking directly for
member functions or free functions in your namespace, it always checks for a valid definition of seqan3::to_rank which
in turn resolves to different possible implementations.

To customise one of our customisation points, do one of the following:

  * Provide functionality as members *or* inside the namespace of your type (will be picked up via
    [argument dependent lookup][https://en.cppreference.com/w/cpp/language/adl).
  * If you adapt a third party's type and you cannot add to that type's namespace, provide functionality inside the
    namespace `seqan3::custom` (this is the "upload namespace").
  * **Never** add names (types, functions, variables...) to namespace `seqan3` and never explicitly specialise one
    of our templates or overload one of our functions.

More technical background on this topic can be found here:

  1. [Original blog article](http://ericniebler.com/2014/10/21/customization-point-design-in-c11-and-beyond/) with
     background and the first "Niebloid" design.
  2. [A more recent article](https://quuxplusone.github.io/blog/2018/03/19/customization-points-for-functions/) that
     also introduces the notion of "upload" namespaces.

There is also a proposed [language extension](http://www.open-std.org/jtc1/sc22/wg21/docs/papers/2018/p1292r0.html)
to handle customisation points in a future version of C++.
