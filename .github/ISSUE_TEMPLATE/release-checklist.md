---
name: Release Checklist
about: A checklist used internally for the release process.
title: "#.#.# Release"
labels: ''
assignees: ''

---

- [ ] Update copyright year.
- [ ] Create at least one pre-release s.t. our package maintainers can try out a new release of our library. ([Instructions](#prerelease))
- [ ] Check that the [directory structure](https://github.com/seqan/seqan3/blob/master/doc/setup/quickstart_cmake/index.md) is valid.
- [ ] Check for critical performance regressions (see [here](https://github.com/seqan/seqan3/wiki/cmp_benchmarks:-Example-Usage)).
- [ ] Check [Nightlies](https://cdash.seqan.de/index.php?project=SeqAn3) for critical build failures.
- [ ] Check workarounds in `platform.hpp`. Are they still valid, or can they be limited to specific compiler versions?
- [ ] Check the [Changelog.md](https://github.com/seqan/seqan3/blob/master/CHANGELOG.md) for completeness (including changed headers). ([Instructions](#changelog))
- [ ] Update the index from cppreference.com so that up-to-date documentation links are generated. ([Instructions](#cppreference))

---

- [ ] Freeze branch. ([Instructions](#freeze))
- [ ] Add versioned documentation to docs.seqan.de. ([Instructions](#versioned-docs))
- [ ] Prepare `seqan3-[VERSION]-{Linux,Source}.tar.xz{,.sha265}`. ([Instructions](#packaging))
- [ ] Prepare a release note with notable features, API changes, bugs, and external dependency updates.

---

- [ ] **Tag release on GitHub and attach `seqan3-[VERSION]-{Linux,Source}.tar.xz{,.sha265}` to the release. ([Instructions](#release))**

---

- [ ] Bump version. ([Instructions](#version-bump))
- [ ] Update the version used for [update notifications](https://github.com/OpenMS/usage_plots/blob/master/seqan_versions.txt).
- [ ] Announce release on [Twitter](https://twitter.com/seqanlib).
- [ ] Announce release on [website](https://www.seqan.de).
- [ ] Announce release on mailing list `seqan-dev@lists.fu-berlin.de`.
- [ ] Announce release on [public Gitter channel](https://gitter.im/seqan/Lobby).
- [ ] Notify upstream package maintainers:
  - [ ] [brew](https://github.com/brewsci/homebrew-bio/tree/develop/Formula/seqan%403.rb) (@eseiler)
  - [ ] [macports](https://github.com/macports/macports-ports/tree/master/science/seqan3/Portfile) (@eseiler)
  - [ ] [bio.tools](https://bio.tools/seqan) (@eseiler)
  - [ ] [conda](https://github.com/bioconda/bioconda-recipes/tree/master/recipes/seqan3) (@eseiler)
  - [ ] [conan](https://github.com/conan-io/conan-center-index/tree/master/recipes/seqan3) (@eseiler)
  - [ ] [debian](https://tracker.debian.org/pkg/seqan3) (@mr-c)
  - [ ] [fedora](https://src.fedoraproject.org/rpms/seqan3) (@sagitter)
- [ ] Update release template with the current release tasks.
- [ ] Celebrate :tada: :beer:

---

#### Instructions

<a name="prerelease"></a>
<details><summary>Creating a pre-release</summary><br>

GitHub is not able to create annotated releases (https://github.com/seqan/product_backlog/issues/159), so we have to manually sign the release.
Make sure you have set up [signed commits](https://docs.github.com/en/authentication/managing-commit-signature-verification/signing-commits).
```bash
git checkout release-[VERSION]
git tag -s [VERSION]-rc.[RC] # e.g. 3.1.0-rc.1
git push upstream [VERSION]-rc.[RC]
```

You will need to provide a tag message. Since this is a pre-release, it can be as simple as `Tag 3.1.0-rc.1`.

Now follow the [packaging instructions](#packaging) to create `seqan3-[VERSION]-rc.[RC]-{Linux,Source}.tar.xz{,.sha265}`.

Go to https://github.com/seqan/seqan3/releases and create a new release using the created tag and attach the source packages.

:warning: **Make sure to set the tick for "This is a pre-release"** :warning:

Once again, the release message can be simply something along the lines of:
```
This is the first release candidate for SeqAn 3.0.3

You can find a list of changes in our [changelog](https://docs.seqan.de/seqan/3.0.3/about_changelog.html).
```

Afterwards, bump the succeeding release candidate number in the release branch: [include/seqan3/version.hpp](https://github.com/seqan/seqan3/blob/3.0.2/include/seqan3/version.hpp#L19-L24).

</details>
<a name="cppreference"></a>
<details><summary>Updating cppreference index</summary><br>

Check for [new releases](https://github.com/PeterFeicht/cppreference-doc/releases) and update the link and hash in [test/documentation/seqan3-doxygen.cmake](https://github.com/seqan/seqan3/blob/b0b279689fa65c2431a5162f2d8acc3ca663f72d/test/documentation/seqan3-doxygen.cmake#L37).
You can compute the hash via `wget -O- <link to html book> | sha256sum`.

</details>
<a name="freeze"></a>
<details><summary>Freezing the release branch</summary><br>

- Make sure all PRs that should be merged are merged.
- Set `SEQAN3_RELEASE_CANDIDATE` to `0` [include/seqan3/version.hpp](https://github.com/seqan/seqan3/blob/3.0.2/include/seqan3/version.hpp#L19-L24).
- This should be the last commit before the release.

</details>
<a name="versioned-docs"></a>
<details><summary>Creating versioned documentation</summary><br>

1. Checkout the release tag and build documentation.
2. Create a #.#.# directory for the release in `/web/docs.seqan.de/htdocs/seqan/`
3. Copy everything from the build (`doc_usr/html/*`) into the directory.
4. Alter the file `/web/docs.seqan.de/htdocs/seqan3.html` with a link to the new documentation build.

</details>
<a name="packaging"></a>
<details><summary>Creating source packages</summary><br>

Use a new clone of the repository.
```bash
git clone https://github.com/seqan/seqan3.git
cd seqan3
git checkout release-[VERSION] # version/branch to pack
git submodule update --init

mkdir ../package-build
cd ../package-build

cmake ../seqan3 # configure
cpack # builds binary package, e.g. seqan3-[VERSION]-Linux.tar.xz{,.sha265}
cmake --build . --target package_source # builds source package, e.g. seqan3-[VERSION]-Source.tar.xz{,.sha265}
```

Note: Do not use `git clone --recurse-submodules https://github.com/seqan/seqan3.git` because it will recursively pull sub-submodules!

</details>
<a name="changelog"></a>
<details><summary>Checking the changelog</summary><br>

- List all supported compiler, also add to https://docs.seqan.de/seqan/3-master-user/about_api.html#autotoc_md35.
- Check that all links are consistent, e.g., `[\#2540](https://github.com/seqan/seqan3/pull/2538)`:
  - Search `(\[\\#)(\d+)(\]\(.+?)(\d+)(\))` and replace `$1$2$3$2$5` (i.e., replace link issue-id by the displayed id).

</details>
<a name="release"></a>
<details><summary>Creating a release</summary><br>

GitHub is not able to create annotated releases (https://github.com/seqan/product_backlog/issues/159), so we have to manually sign the release.
Make sure you have set up [signed commits](https://docs.github.com/en/authentication/managing-commit-signature-verification/signing-commits).
```bash
git checkout release-[VERSION]
git tag -s [VERSION]
git push upstream [VERSION]
```

You will need to provide a tag message. We use the first sentences of the release note:

E.g. (see https://github.com/seqan/seqan3/tags)
```
SeqAn 3.0.2 Release


Despite all circumstances, we are excited to present a new update of our SeqAn library.
We present some great new features and also a lot of usability improvements.
Among others, this release will fully comply with the final C++-20 standard.

:warning: In this release we harmonised the algorithm configurations for a better user experience.
This, much like 2020, will break a lot of code. But rest assured that the changes are easy to apply and are worth every bit. :smile:

You can find a comprehensive list of the changes in our [changelog](https://docs.seqan.de/seqan/3.0.2/about_changelog.html).
```

</details>
<a name="version-bump"></a>
<details><summary>Bumping the version</summary><br>

- Bump succeeding version number in the master branch: [include/seqan3/version.hpp](https://github.com/seqan/seqan3/blob/3.0.2/include/seqan3/version.hpp#L19-L24).
- The `SEQAN3_RELEASE_CANDIDATE` must be set to `1` as `0` indicates a stable release.
- Bump the latest stable version number of the API-Stability test in master: [test/api_stability/CMakeLists.txt](https://github.com/seqan/seqan3/blob/3.0.3/test/api_stability/CMakeLists.txt#L10).

</details>
