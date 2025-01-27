<!-- SPDX-FileCopyrightText: 2006-2025 Knut Reinert & Freie Universität Berlin
     SPDX-FileCopyrightText: 2016-2025 Knut Reinert & MPI für molekulare Genetik
     SPDX-License-Identifier: CC-BY-4.0
-->

# SeqAn3 API Stability

This test suite will check whether the current version of seqan3 (e.g. git checkout) can be used to build the test cases
of the latest stable version of seqan3 (https://github.com/seqan/seqan3/releases/latest).

By seeing which of our own previous test cases break, we can estimate how difficult it is to upgrade an external app and
adapt our changelog accordingly.

### Usage

```bash
# Note that the build directory must not be a subdirectory of a git repository (a directory containing .git)
# Otherwise, the patches cannot be applied correctly
cd <build_dir> # e.g. seqan3-builds/api_stability/

# we recommend to disable deprecation warnings in CI
cmake -DCMAKE_CXX_FLAGS="-DSEQAN3_DISABLE_DEPRECATED_WARNINGS=1" <seqan3_git_checkout>

# if you want to build with 40 threads, you need to use CMAKE_BUILD_PARALLEL_LEVEL to specify the threads
# a normal `make -j 40` or `cmake --build . -j 40` will not work
CMAKE_BUILD_PARALLEL_LEVEL=40 cmake --build .
```

### How does this work?

In the following, we will use `3.3.0` as the latest stable release (LSR).

The workflow entails these steps:
1. Download and unzip the LSR.
2. Remove everything from the LSR except the tests in `<seqan3_lsr_source>/test/`.
3. Apply patches from the current git version on the LSR, i.e., patches found in
   `<seqan3_git_source>/test/api_stability/3.3.0`.
  * These patches will only apply changes on the tests.
  * This step is necessary as some of our tests also test non-public / non-stable API.
3. Use the current version of `find_package (SeqAn3)` found in `<seqan3_git_source>/cmake`.
4. Configure the tests `<seqan3_lsr_source>/test/unit` with our current seqan3 header-library
   (e.g. `<seqan3_git_source>/include` and `<seqan3_git_source>/submodules/*/include`).
  * But use the test cases from the LSR.
5. Build all tests.
6. Run all tests.
7. If any of the above steps fail, the complete build will fail.

We can apply this workflow to every kind of test (e.g., `test/unit`, `test/snippet`, `test/performance`) that we have.
Currently, we only test `test/unit` and `test/snippet` as both of them test most of our public API.

### When to create patches and what is the naming scheme?

Patches should only be created if API changed in a way that we cannot give an "automatic" upgrade path by, for example,
defining a type alias, a compatible function that perfectly forwards to the new API or some other means that provide the
old API which then uses the new API. Please try your hardest to provide "compatible" API, as this is the easiest for the
user.

There are two categories of API changes that are reflected in the patches (prefixes in the commit message):

* `[API]` are patches that changed the public API in some way (regardless of whether this is `stable` or `experimental`
  API), like the behaviour, the return type or the call convention. These patches should always have a `CHANGELOG.md`
  entry explaining what the change did. `[API]` changes are further divided into these subcategories:
  * `[API] [BREAKAGE]` An API that existed and has a completely different meaning now, such a change MUST have a really
    good reason and should be avoided or documented extremely well.
  * `[API] [FIX]` The change was necessary, as it fixed a defect. Some examples:
    * A function does not throw anymore / will throw now.
    * Format of output (e.g., `std::cout`) changed.
    * CTAD changed (e.g., `some_class{0}` will now be of type `some_class<int>` instead of `some_class<double>`).
    * The return type changed / order within template parameters changed.
  * `[API] [WRONG-USAGE]` An API was misused in the test or tested some invalid scenarios.
* `[NOAPI]` are patches that change internal (e.g., `seqan3::detail` or `private` member) or non-public API (e.g.,
  marked `\noapi`). These changes do not need an upgrade path. We differentiate between the following subcategories:
  * `[NOAPI] [INCLUDE]` Some header file was moved, or an unnecessary include was removed in some header,
    which caused a test to explicitly include that header. As there are tools like [Include What You
    Use](https://github.com/include-what-you-use/include-what-you-use), we do not consider these changes as API
    stability issues. We will still try to add an upgrade path if possible (like when we split `/include/seqan3/core`
    and `/include/seqan3/utility`) and add `CHANGELOG.md` entries for every header change.
  * `[NOAPI] [DETAIL]` Some `seqan3::detail` code changed, like a rename, or template parameter change.
  * `[NOAPI] [DEPRECATED]` Some deprecated code was removed.
  * `[NOAPI] [BREAKAGE]` The same API has a considerable different behaviour. E.g., it accepted some
    input before, but now it would segfault on the same input.


### How to create patches?

In the following, we will use `3.3.0` as the latest stable release (LSR).

Create a new branch based on the LSR and apply all existing patches in `<seqan3_git_source>/test/api_stability/3.3.0`.

```
cd <seqan3_git_source>

# assume that your current branch you are working on is test/api
git checkout test/api

# copy over patches to a tmp directory (`git am` seems to not support applying patches onto a different branch)
mkdir -p /tmp/seqan3-api-stability-patches
cp test/api_stability/3.3.0/*.patch /tmp/seqan3-api-stability-patches

# create a new branch based on the LSR and switch to it
git checkout -b api-stability-patches 3.3.0

# apply all patches onto 3.3.0 (--keep-non-patch will keep `[NOAPI]` tags in the commit message)
git am --keep-non-patch /tmp/seqan3-api-stability-patches/*.patch

# clean up applied patches
rm /tmp/seqan3-api-stability-patches/*.patch
```

Now re-apply the commit(s) that changed the API.

```
git cherry-pick <commit-SHA>
# To apply a range of commits, e.g. from a PR, you can use `git chery-pick <first-commit-SHA>..<second-commit-SHA>`

# fix possible merge conflicts
git cherry-pick --continue

# remove everything except test/unit, test/snippet, and doc/**.cpp

# commit everything to the current cherry-picked commit
git commit --amend
```

Alternatively, you can also first resolve the merge conflicts and use `git reset` to reset the index to the state
before cherry-picking. Then, you can remove unecessary files and commit with an appropriate commit message.

You will most likely get merge conflicts. You can ignore all merge conflicts in any folder other than
`<seqan3_lsr_source>/test/unit`, `<seqan3_lsr_source>/test/snippet`, and `<seqan3_lsr_source>/doc`.

It is also important that you double-check whether the patch only contains changes that apply to
`<seqan3_lsr_source>/test/{unit, snippet}` and `<seqan3_lsr_source>/doc/**.cpp`, any other change must be discarded.
New files for tests must not be added.

After that, we can export all patches.

```
# export all patches since 3.3.0
git format-patch 3.3.0

# move them to tmp directory
mv *.patch /tmp/seqan3-api-stability-patches
```

Now change to your branch that you were working on and check-in the patches.

```
git checkout test/api

cp /tmp/seqan3-api-stability-patches/*.patch test/api_stability/3.3.0/

rm -fdr /tmp/seqan3-api-stability-patches
git branch -D api-stability-patches

# add new patches
git add test/api_stability/3.3.0/

# commit changes
git commit
```

Before pushing, try whether the patches work as intended.
