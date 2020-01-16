---
name: Release Checklist
about: A checklist used internally for the release process.
title: "#.#.# Release"
labels: ''
assignees: ''

---

- [ ] Check [Nightlies](http://cdash.seqan.de/index.php?project=SeqAn3) for critical build failures
- [ ] Check Changelog.md for completeness
- [ ] Prepare `seqan-[VERSION]-with-submodules.tar.gz`
- [ ] Tag release on Github and attach `seqan-[VERSION]-with-submodules.tar.gz` to the release
- [ ] Bump succeeding version number
- [ ] update the index from cppreferences.com so that up-to-date documentation links are generated.
- [ ] add versioned documentation build on docs.seqan.de
- [ ] Announce release on twitter
- [ ] Announce release on www.seqan.de
- [ ] Announce release on mailing list `seqan-dev@lists.fu-berlin.de`
- [ ] Announce release on public Gitter channel [![Gitter](https://badges.gitter.im/seqan/Lobby.svg)](https://gitter.im/seqan/Lobby?utm_source=badge&utm_medium=badge&utm_campaign=pr-badge)
- [ ] Celebrate :tada: :beer:
