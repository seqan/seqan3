# Stability and long-term promises {#about_api}

<!-- SPDX-FileCopyrightText: 2006-2025 Knut Reinert & Freie Universität Berlin
     SPDX-FileCopyrightText: 2016-2025 Knut Reinert & MPI für molekulare Genetik
     SPDX-License-Identifier: CC-BY-4.0
-->

[TOC]

SeqAn3 adheres to [semantic versioning](https://semver.org) and provides a stable API within
one major release (all of `seqan-3.*.*`) unless otherwise noted.

There is no ABI stability.

Many of the rules on this page are derived from the [abseil library](https://abseil.io/about/compatibility).

# API stability {#api_stability}

In general, you can expect all **stable** entities within the `seqan3` namespace to be usable throughout the entire
release cycle of SeqAn3, i.e. if you write an application that includes `seqan-3.1.1`, you should be able
to compile it against `seqan-3.4.5` without errors.

Exceptions to the previous rule:
  -# **Do not depend upon internal details.** If something is in a namespace or filename/path that includes the word
     `detail` or if it is private to a class type, you are not allowed to depend upon it in any way. More generally,
     if it is not part of the *user documentation*, it is not part of the API.
  -# **Do not depend on entities marked as "experimental" or "no-api"**. Major changes to the library like new modules
     are often first marked as experimental within a minor release cycle. This means we do not guarantee
     stability until the next minor release happens and the experimental flag is removed in favour of the stable flag.
     Some entities in namespace `seqan3` are permanently marked as "no-api" which designates them as subject to
     unannounced change. In addition, all entities in **core** and **utility** are always marked "no-api". This is the
     case for auxiliary data structures (usually in `seqan3::detail`) that are needed for the public API documentation
     or because they are considered a nice-to-have feature for users.
  -# **Do not depend on the *signatures* of SeqAn APIs.** In particular, you may not take the address of APIs in SeqAn
     and you may not use metaprogramming tricks to depend on those signatures. We reserve the right to:
     * Add new names to namespace `seqan3` and any sub-namespaces
     * Add new member functions to types in namespace `seqan3`
     * Add new overloads to existing functions
     * Add new default arguments to functions and templates
     * Change return-types of functions in compatible ways (void to anything, numeric types in a widening fashion, etc.).
  -# **Do not forward declare SeqAn APIs.** This is actually a sub-point of the previous one, but can be
     surprising. Any refactoring that changes template parameters, default parameters, or namespaces will be a breaking
     change in the face of forward-declarations.
  -# **Avoid unnecessary dependency on Argument-Dependent Lookup (ADL) when calling SeqAn APIs.** Some APIs are designed
     to work via ADL (e.g. `operator<<` for iostreams, unqualified `swap` in generic code, etc.) For most APIs, however,
     ADL is not part of the design. Calling functions from namespace `seqan3` via ADL, unless that is explicitly
     intended as part of the design, should be avoided.
  -# **Include What You Use.** We may make changes to the internal include-graph for SeqAn headers - if you use an
     API, please include the relevant header file directly.
  -# **Do not make unqualified calls in the global namespace.** A call like `f(a);` for a function `f` in the global
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
As a special case of point 2, **the entire 3.0.* release is not stable**.
We will do our best not to break things, but similar to the releases of GCC, we start labelling entities as stable
starting from `seqan-3.1.0`. Please refer to the documentation of the individual entities for their stability status.

# ABI stability

**No guarantees of any kind** are provided in regard to the in-memory representation of our data structures.
We intend for SeqAn to be built from source. The internal layout of our types may change at any point, without notice.
In particular, building SeqAn in the presence of different C++ standard library types may change SeqAn types,
especially for pre-adopted types in the `std` module — these will resolve to standard library types as of C++20.

We do promise, however, that SeqAn3 data structures¹ serialised to disk are de-serialisable by later versions of SeqAn3.
The reverse is not true, however, and re-serialising a data structure serialised by a previous release might result in
a different (updated) on-disk format.

<small>¹ Iff they are API stable.</small>

# Platform stability {#platform_stability}

The main requirement for SeqAn3 is that your operating system provides one of the compilers supported by us.
In general, we only support the latest three major compiler versions.
We currently support the following compilers on 64-bit operating systems with little-endian CPU architectures:
  * GCC11, GCC12, GCC13

\note Only the most recent minor release of a compiler is guaranteed to be supported, e.g. when `gcc-11.4` is released,
we may drop support for `gcc-11.3`. Since all platforms with an older version receive minor release updates,
this should not be a problem.

We promise to support the above compilers in the latest release of SeqAn3, or until all the following
operating systems provide a newer supported compiler:

| Operating System             | Supported Releases¹                    |
|------------------------------|----------------------------------------|
| RedHat Enterprise Linux      | the latest release ²                   |
| CentOS Linux                 | the latest release ²                   |
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

More platforms and compilers will be added during the SeqAn3 lifecycle, but some will be marked as
*experimentally supported* which means support may be dropped again later on.

**We promise to provide good forward-compatibility with the C++ standard.** And we will strive to fix any warnings that
are added by newer versions of a supported compiler.

<small>¹ [This site](https://linuxlifecycle.com) provides a good overview of what the current release and its
lifecycle is.</small><br>
<small>² We consider CentOS 7 / RedHat Enterprise Linux (RHEL) 7 as being community-supported. That means issues and
patches are welcome, but we do not actively test for those operating systems. See this related
[issue](https://github.com/seqan/seqan3/issues/2244).</small>

# Dependencies

We may add additional dependencies in future minor releases of SeqAn or raise the required versions of current
dependencies.
However, we promise not to add *required* dependencies that result in new licensing obligations or link-time
dependencies.
