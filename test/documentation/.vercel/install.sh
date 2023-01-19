#!/usr/bin/env bash
set -exo pipefail

DOXYGEN_VERSION=1.9.6
SOURCE_DIR=`pwd`
CACHE_DIR="${SOURCE_DIR}/node_modules"

# Install dependencies.
yum --assumeyes --quiet install wget cmake3 flex bison xz graphviz

# Download doxygen.
mkdir -p ${CACHE_DIR}/doxygen-download
wget --retry-connrefused --waitretry=30 --read-timeout=30 --timeout=30 --tries=20 --no-clobber --quiet --directory-prefix=${CACHE_DIR}/doxygen-download/ https://github.com/doxygen/doxygen/releases/download/Release_${DOXYGEN_VERSION//./_}/doxygen-${DOXYGEN_VERSION}.src.tar.gz
tar -C ${CACHE_DIR}/doxygen-download -zxf ${CACHE_DIR}/doxygen-download/doxygen-${DOXYGEN_VERSION}.src.tar.gz

# Configure doxygen.
cd ${CACHE_DIR}/doxygen-download/doxygen-${DOXYGEN_VERSION}
mkdir -p build && cd build
if ! cmake3 -G "Unix Makefiles" .. 1>/dev/null; then
    rm -rf *
    cmake3 -G "Unix Makefiles" .. 1>/dev/null
fi

# Build doxygen.
make -j 4 1>/dev/null

# Install doxygen
make install DESTDIR=${CACHE_DIR}/doxygen-${DOXYGEN_VERSION} 1>/dev/null

# Symlink doxygen.
ln -s ${CACHE_DIR}/doxygen-${DOXYGEN_VERSION}/usr/local/bin ${SOURCE_DIR}/doxygen-bin
