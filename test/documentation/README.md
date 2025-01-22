<!-- SPDX-FileCopyrightText: 2006-2025 Knut Reinert & Freie Universität Berlin
     SPDX-FileCopyrightText: 2016-2025 Knut Reinert & MPI für molekulare Genetik
     SPDX-License-Identifier: CC-BY-4.0
-->

# seqan3 documentation

Currently, we can only build the documentation on *nix systems.

We offer two versions of our documentation, one intended for the user (doc_usr) and one intended for the
library-developer and maintainer (doc_dev) of seqan3, which contains the documentation of internals.

## How to configure:

```bash
mkdir <seqan3-build-dir>/documentation
cd <seqan3-build-dir>/documentation

cmake <seqan3-dir>/test/documentation
```

## How to build:

Prerequisites: configuring the documentation

```bash
cd <seqan3-build-dir>/documentation

# build user and developer documentation
cmake --build .

# or build user documentation only
cmake --build . --target doc_usr

# or build developer documentation only
cmake --build . --target doc_dev
```

## How to test:

Prerequisites: building the documentation

```bash
cd <seqan3-build-dir>/documentation

# test user and developer documentation
ctest --output-on-failure --progress

# or test user documentation only
ctest -R doc_usr --output-on-failure --progress

# or test developer documentation only
ctest -R doc_dev --output-on-failure --progress
```

## How to install:

Prerequisites: building the documentation

Our installation uses GNU standard installation directories provided by cmake
[GNUInstallDirs](https://cmake.org/cmake/help/v3.19/module/GNUInstallDirs.html#module:GNUInstallDirs).

That means the html documentation will be installed to `<DESTDIR>/<INSTALL_PREFIX>/<INSTALL_DOCDIR>/html/` where 
`<INSTALL_PREFIX>` is typically `usr/local` and `<INSTALL_DOCDIR>` is `share/doc/seqan3`.

That means our html documentation will be installed to `<DESTDIR>/usr/local/share/doc/seqan3/html`.

You can override the paths by specifying `-DCMAKE_INSTALL_PREFIX` and `-DCMAKE_INSTALL_DOCDIR` during configuration time

```bash
# Note that -DCMAKE_INSTALL_DOCDIR="." is important, otherwise it will install it to `share/doc/seqan3`
cmake -DCMAKE_INSTALL_PREFIX="" -DCMAKE_INSTALL_DOCDIR="." <seqan3-dir>/test/documentation
```

which will install the documentation to `<DESTDIR>/`, i.e. without any prefixes.

```bash
cd <seqan3-build-dir>/documentation

# per default we only install the user documentation
# --prefix will set <INSTALL_PREFIX> to export
# Note that combining DESTDIR with --prefix might result in weird output
#
# ./export/<INSTALL_DOCDIR>/html
cmake --install . --prefix export # (since cmake 3.15) or

# ./export/<INSTALL_PREFIX>/<INSTALL_DOCDIR>/html/
DESTDIR="export" cmake --install . # (since cmake 3.15)

# the user documentation can be installed via
#
# ./export/<INSTALL_DOCDIR>/html
cmake --install . --prefix export --component doc # (since cmake 3.15) or

# ./export/<INSTALL_PREFIX>/<INSTALL_DOCDIR>/html
DESTDIR="export" cmake -DCOMPONENT=doc -P cmake_install.cmake # (for older cmake versions)

# the developer documentation can be installed via
# --prefix will set <INSTALL_PREFIX> to export-dev
#
# ./export-dev/<INSTALL_DOCDIR>/html
cmake --install . --prefix export-dev --component doc-dev # (since cmake 3.15) or

# ./export-dev/<INSTALL_PREFIX>/<INSTALL_DOCDIR>/html
DESTDIR="export-dev" cmake -DCOMPONENT=doc-dev -P cmake_install.cmake # (for older cmake versions)
```

## How to package:

Prerequisites: building the documentation

```bash
cd <seqan3-build-dir>/documentation

# Create user documentation package, i.e. seqan3-3.0.2-Linux-doc.tar.gz{,.sha256}
# The package will contain the following structure <INSTALL_PREFIX>/<INSTALL_DOCDIR>/html
cpack

# Create user and developer documentation package, i.e.
# seqan3-3.0.2-Linux-{doc,doc-dev}.tar.gz{,.sha256}
cpack -D CPACK_COMPONENTS_ALL="doc;doc-dev"
```
