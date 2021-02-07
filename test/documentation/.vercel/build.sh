#!/usr/bin/env bash
set -x
set -e

SOURCE_DIR=`pwd`
BUILD_DIR=`pwd`/build
EXPORT_DIR=`pwd`/export

PATH="$PATH:${SOURCE_DIR}/doxygen-bin"

echo $PATH
ls -al ${SOURCE_DIR}/doxygen-bin
${SOURCE_DIR}/doxygen-bin/doxygen --version

cmake3 --version
doxygen --version

mkdir -p "${BUILD_DIR}"
cd "${BUILD_DIR}"

cmake3 -DCMAKE_INSTALL_PREFIX="" -DCMAKE_INSTALL_DOCDIR="." "${SOURCE_DIR}/test/documentation"

# cmake3 --build . --target download-cppreference-doxygen-web-tag
cmake3 --build .

mkdir -p "${EXPORT_DIR}/usr/" "${EXPORT_DIR}/dev/"

DESTDIR="${EXPORT_DIR}/usr/" cmake3 -DCOMPONENT=doc -P cmake_install.cmake
DESTDIR="${EXPORT_DIR}/dev/" cmake3 -DCOMPONENT=doc-dev -P cmake_install.cmake

cp "${SOURCE_DIR}/test/documentation/.vercel/index.html" "${EXPORT_DIR}/index.html"
