#!/usr/bin/env bash
set -x
set -e

DOXYGEN_VERSION=1.9.1

SOURCE_DIR=`pwd`
CACHE_DIR="${SOURCE_DIR}/node_modules"

yum -y install graphviz texlive texlive-utils wget cmake3 flex bison xz # doxygen

mkdir -p ${CACHE_DIR}/doxygen-download
wget --retry-connrefused --waitretry=30 --read-timeout=30 --timeout=30 --tries=20 --no-clobber --quiet --directory-prefix=${CACHE_DIR}/doxygen-download/ https://sourceforge.net/projects/doxygen/files/rel-${DOXYGEN_VERSION}/doxygen-${DOXYGEN_VERSION}.src.tar.gz
tar -C ${CACHE_DIR}/doxygen-download -zxf ${CACHE_DIR}/doxygen-download/doxygen-${DOXYGEN_VERSION}.src.tar.gz

pushd ${CACHE_DIR}/doxygen-download/doxygen-${DOXYGEN_VERSION}
mkdir -p build && cd build
if ! cmake3 -G "Unix Makefiles" ..; then
    rm -rf *
    cmake3 -G "Unix Makefiles" ..
fi
make -j 4
make install DESTDIR=${CACHE_DIR}/doxygen-${DOXYGEN_VERSION}
popd
# tar -C ${CACHE_DIR}/ -zxf ${CACHE_DIR}/doxygen-download/doxygen-${DOXYGEN_VERSION}.linux.bin.tar.gz

# Only available in subsequent steps!
ln -s ${CACHE_DIR}/doxygen-${DOXYGEN_VERSION}/usr/local/bin doxygen-bin
