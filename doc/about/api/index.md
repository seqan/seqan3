# Stability and future-proofeness {#about_api}

[TOC]

SeqAn3 adheres to [semantic versioning](https://semvar.org) and provides stable interface within
one major release (all of `seqan-3.*.*`) unless otherwise noted.
Being a header-only library, it makes no promises about ABI stability, though.
Fore more details, see below.

# API stability

In general, you can expect all entities within the `seqan3` namespace to be usable throughout the entire
release cycle of SeqAn3, i.e. if you write an application that includes `seqan-3.1.1`, you should be able
to re-compile it against `seqan-3.4.5` without errors.

Notable exceptions to this rule:
  1. There are no guarantees about entities in the `seqan3::detail` namespace which is why this is excluded from
     the reguler user documentation.
  2. Major changes to the library like new modules are often marked as experimental within the first minor release
     they appear in. This means we do not guarantee stability until the next minor release happens and the experimental
     flag is removed. (This allows us to gather more user experience with new interfaces and make minor changes if
     necessary.)
  3. Very few entities in namespace `seqan3` are permanently marked as "NOAPI" which designates them as subject to
     unannounced change. This is usually the case for auxiliary data structures that would be part of
     `seqan3::detail` were they not needed to generate correct API documentation for entities in `seqan3`.

All of the above cases will be distinctly visible in the documentation.

\warning
As a special case of point 2. **the entire 3.0.* release is marked as experimental**.
We will do our best not to break things, but similar to the releases of GCC, you should consider `seqan-3.1.0` as
the first API-stable release of SeqAn3.

\note
We only promise that your applications build *error-free* against future releases of SeqAn3, in certain
situations you might receive *warnings*, though.
In particular we will mark certain entities as `[[deprecated]]` when a superior interface is added.
Being API-stable means you *can* ignore these warnings, but we do recommend following update instructions.

# ABI stability

**No guarantees of any kind** are provided in regard to the in-memory representation of our data structures.
We do promise, however, that SeqAn3 data structures¹ serialised to disk are de-serialisable by later versions of SeqAn3.
The reverse is not true, however, and re-serialising a data structure serialised by a previous release might result in
a different (updated) on-disk format.

<small>¹ Iff they are API stable.</small>

# Platform stability

The following operating systems/compilers are guaranteed to be supported until the last release of SeqAn3:
  * GNU/Linux, macOS: GCC7, GCC8, GCC9

\note
Only the most recent minor release of a compiler is guaranteed to be supported, e.g. when `gcc-7.5` is released,
we may drop support for `gcc-7.4`.
Since all platforms with an older version receive minor release updates, this should not be a problem.

More platforms and compilers will be added during the SeqAn3 lifecycle, but some will be marked as
*experimentally supported* which means support may be dropped later on.

# Dependencies

We may add additional dependencies in future minor releases of SeqAn or raise the required versions of current
dependencies.
However, we promise not to add *required* dependencies that result in new licensing obligations or link-time
dependencies.
