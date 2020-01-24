---
name: Release Checklist
about: A checklist used internally for the release process.
title: "#.#.# Release"
labels: ''
assignees: ''

---

- [ ] Check for critical performance regressions.
- [ ] Check [Nightlies](http://cdash.seqan.de/index.php?project=SeqAn3) for critical build failures.
- [ ] Check the [Changelog.md](https://github.com/seqan/seqan3/blob/master/CHANGELOG.md) for completeness.
- [ ] Update the index from cppreference.com so that up-to-date documentation links are generated.
  <details><summary>Click for detailed steps</summary><br>

  Check for [new releases](https://github.com/PeterFeicht/cppreference-doc/releases) and update the link and hash
  in the [code base](https://github.com/seqan/seqan3/blob/b0b279689fa65c2431a5162f2d8acc3ca663f72d/test/documentation/seqan3-doxygen.cmake#L37).

  For the hash do

  ```
  wget -O- <link to html book> | sha256sum
  ```

  </details>
- [ ] Add versioned documentation build on docs.seqan.de.
  <details><summary>Click for detailed steps</summary>

  1. Build the documentation locally

  2. Create a #.#.# directory for the release in `/web/docs.seqan.de/htdocs/seqan/`

  3. Copy everything from the local build (`doc_usr/html/*`) into the directory.

  4. Alter the file `/web/docs.seqan.de/htdocs/index.html` with a link to the new documentation build.

  </details>
- [ ] Prepare `seqan-[VERSION]-with-submodules.tar.gz`.
  <details><summary>Click for detailed steps</summary><br>

  Make a fresh clone of the repository and recursively delete the `.git` directories.
  ```
  git clone https://github.com/seqan/seqan3.git
  cd seqan3
  git submodule init
  git submodule update
  rm -rdf .git submodules/sdsl-lite/.git submodules/range-v3/.git submodules/cereal/.git
  cd ..
  tar -czf seqan-<#.#.#>-with-submodules.tar.gz seqan3
  ```

  Note: Do not do `git clone -recurse-submodules https://github.com/seqan/seqan3.git` because it will recursively
  pull sub-submodules which we do not want!

  </details>
- [ ] Prepare a release note with notable features, API changes, bugs, and external dependency updates.

---

- [ ] **Tag release on Github and attach `seqan-[VERSION]-with-submodules.tar.gz` to the release.**

---

- [ ] Bump succeeding version number. ([include/seqan3/version.hpp](https://github.com/seqan/seqan3/blob/6391aaa014db4d04b8220a2041014295c4134f33/include/seqan3/version.hpp#L19-L24))
- [ ] Update the SeqAn version [here](https://github.com/OpenMS/usage_plots/blob/master/seqan_versions.txt) to ensure
      that the server in Tübingen, which takes care of the argument parser update notifications, is aware of it
      (@smehringer).
- [ ] Announce release on twitter.
- [ ] Announce release on www.seqan.de.
- [ ] Announce release on mailing list `seqan-dev@lists.fu-berlin.de`.
- [ ] Announce release on public Gitter channel [![Gitter](https://badges.gitter.im/seqan/Lobby.svg)](https://gitter.im/seqan/Lobby?utm_source=badge&utm_medium=badge&utm_campaign=pr-badge).
- [ ] Notify upstream package maintainers
  - brew (@rrahn)
  - macports (@rrahn)
  - conda (@eseiler)
  - debian (@mr-c)
- [ ] Celebrate :tada: :beer:
