# Stability and future-proofness {#about_api}

[TOC]

Summary: SeqAn3 adheres to [semantic versioning](https://semvar.org) and provides a stable API within
one major release (all of `seqan-3.*.*`) unless otherwise noted.
There is no ABI stability.

Many of the rules on this page are derived from the [abseil library](https://abseil.io/about/compatibility).

# API stability

In general, you can expect all entities within the `seqan3` namespace to be usable throughout the entire
release cycle of SeqAn3, i.e. if you write an application that includes `seqan-3.1.1`, you should be able
to re-compile it against `seqan-3.4.5` without errors.

Exceptions to the previous rule:
  1. **Do not depend upon internal details.** If something is in a namespace or filename/path that includes the word
     `detail` or if it is private to a class type, you are not allowed to depend upon it any way. More generally,
     if it is not part of the *user documentation*, it is not part of the API.
  2. **Do not depend on entities marked as "experimental" or "NOAPI"**. Major changes to the library like new modules
     are often marked as experimental within the first minor release they appear in. This means we do not guarantee
     stability until the next minor release happens and the experimental flag is removed. Very few entities in namespace
     `seqan3` are permanently marked as "NOAPI" which designates them as subject to unannounced change. This is usually
     the case for auxiliary data structures that would be part of `seqan3::detail` were they not needed to generate
     correct API documentation for entities in `seqan3`.
  4. **Do not depend on the *signatures* of SeqAn APIs.** In particular you may not take the address of APIs in SeqAn
     and you may not use metaprogramming tricks to depend on those signatures. We reserve the right to:
     * Add new names to namespace `seqan3` and any sub-namespaces
     * Add new member functions to types in namespace `seqan3`
     * Add new overloads to existing functions
     * Add new default arguments to functions and templates
     * Change return-types of functions in compatible ways (void to anything, numeric types in a widening fashion, etc.).
  5. **Do not forward declare SeqAn APIs.** This is actually a sub-point of the previous one, but can be
     surprising. Any refactoring that changes template parameters, default parameters, or namespaces will be a breaking
     change in the face of forward-declarations.
  6. **Avoid unnecessary dependency on Argument-Dependent Lookup (ADL) when calling SeqAn APIs.** Some APIs are designed
     to work via ADL (e.g. `operator<<` for iostreams, unqualified `swap` in generic code, etc.) For most APIs, however,
     ADL is not part of the design. Calling functions from namespace `seqan3` via ADL, unless that is explicitly
     intended as part of the design, should be avoided.
  7. **Include What You Use.** We may make changes to the internal include-graph for SeqAn headers - if you use an
     API, please include the relevant header file directly.
  8. **Do not make unqualified calls in the global namespace.** A call like `f(a);` for a function `f` in the global
     namespace can become ambiguous if/when we add `seqan3::f` (especially if `a` is a SeqAn type). We generally do
     not recommend you use the global namespace for anything. If you must, please qualify any call that accepts a type
     provided by SeqAn.

In the very rare case that we feel a breaking change to the API is unavoidable, we promise to provide the following
update path:
  * A new minor version is introduced that supports the new API.
  * The old API is still supported, but might be marked `[[deprecated]]`.
  * We provide very extensive documentation on the change or an automated tool that converts your code from using the
    old API to using the new API.
  * The next (or another future) minor version removes the old API.

\warning
As a special case of point 2. **the entire 3.0.* release is marked as experimental**.
We will do our best not to break things, but similar to the releases of GCC, you should consider `seqan-3.1.0` as
the first API-stable release of SeqAn3.

# ABI stability

**No guarantees of any kind** are provided in regard to the in-memory representation of our data structures.
We intend for SeqAn to be built from source. The internal layout of our types may change at any point, without notice.
In particular, building SeqAn in the presence of different C++ standard library types may change SeqAn types,
especially for pre-adopted types in the `std` module — these will resolve to standard library types as of C++20.

We do promise, however, that SeqAn3 data structures¹ serialised to disk are de-serialisable by later versions of SeqAn3.
The reverse is not true, however, and re-serialising a data structure serialised by a previous release might result in
a different (updated) on-disk format.

<small>¹ Iff they are API stable.</small>

# Platform stability

The main requirement for being able to use SeqAn is a supported compiler and in general all operating systems that
provide a supported compiler should be supported. We currently support the
following compilers on 64-bit operating systems with little-endian CPU
architectures only:
  * GCC7, GCC8, GCC9

We promise to support all of the above compilers until the last release of SeqAn3, or until all of the following
operating systems provide a newer supported compiler:

| Operating System             | Supported Releases¹                    |
|------------------------------|----------------------------------------|
| RedHat Enterprise Linux      | the latest release                     |
| CentOS Linux                 | the latest release                     |
| SUSE Linux Enterprise Server | the latest release                     |
| Debian GNU/Linux             | "stable" and "old-stable"              |
| Ubuntu Linux                 | the two latest LTS releases            |
|                              |                                        |
| macOS                        | the two latest releases                |
|                              |                                        |
| FreeBSD                      | the latest stable release              |

"Support" in this context does not imply that we package SeqAn for these operating systems (although we do in some
cases), it merely states that you should be able to build applications that depend on SeqAn on the given platforms.
This implies the availability of a supported compiler in the default package repositories or via easy-to-use
third party services.

\note
Only the most recent minor release of a compiler is guaranteed to be supported, e.g. when `gcc-7.5` is released,
we may drop support for `gcc-7.4`.
Since all platforms with an older version receive minor release updates, this should not be a problem.

More platforms and compilers will be added during the SeqAn3 lifecycle, but some will be marked as
*experimentally supported* which means support may be dropped again later on.

**We promise to provide good forward-compatibility with the C++ standard.** And we will strive to fix any warnings that
are added by newer versions of a supported compiler.

<small>¹ [This site](https://linuxlifecycle.com) provides a good overview of what the current release and its
lifecycle is.</small>

# Dependencies

We may add additional dependencies in future minor releases of SeqAn or raise the required versions of current
dependencies.
However, we promise not to add *required* dependencies that result in new licensing obligations or link-time
dependencies.
