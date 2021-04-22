#!/usr/bin/env bash
set -x
set -e

sudo apt-get install texlive-font-utils ghostscript texlive-latex-extra graphviz libclang-9-dev libclang-cpp9 # graphviz for dot, latex to parse formulas, libclang for doxygen
mkdir -p /tmp/doxygen-download
wget --retry-connrefused --waitretry=30 --read-timeout=30 --timeout=30 --tries=20 --no-clobber --quiet --directory-prefix=/tmp/doxygen-download/ https://sourceforge.net/projects/doxygen/files/rel-${DOXYGEN_VERSION}/doxygen-${DOXYGEN_VERSION}.linux.bin.tar.gz
tar -C /tmp/ -zxf /tmp/doxygen-download/doxygen-${DOXYGEN_VERSION}.linux.bin.tar.gz
echo "/tmp/doxygen-${DOXYGEN_VERSION}/bin" >> $GITHUB_PATH # Only available in subsequent steps!
