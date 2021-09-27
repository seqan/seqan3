#!/usr/bin/env bash

case $# in
    1 )
        package_name=$1
        brew_package_name=$1 ;;
    2 )
        package_name=$1
        brew_package_name="${1}@${2}" ;;
    * )
        echo "Usage: install_via_brew.sh <package_name> [version]"
        exit 1 ;;
esac

set -e
set -x

# General idea:
# `brew list --versions` returns 1 if package is installed, and 0 otherwise. By combining it with some other command
# via `&&`, we can run the second command depending on whether the package is installed.
# package_name: The package name without a version, e.g. "ccache" and "gcc".
# brew_package_name: The package followed by an optional version, e.g. "ccache" and "gcc@11".
# `--force-bottle` will cause brew to use precompiled binaries. Sometimes brew likes to build compilers from scratch.

# If the package is installed, we unlink the package. brew unlink does not take a version.
brew list --versions $brew_package_name && \
    brew unlink $package_name \
|| true

# If the package is installed, we upgrade. Otherwise, we install it.
brew list --versions $brew_package_name && \
    brew upgrade --force-bottle $brew_package_name \
|| brew install --force-bottle $brew_package_name

# We link the package, i.e. we add symlinks into /usr/local/bin. They requested version of the package should superseed
# other installed versions (`--overwrite`). For example, when requesting gcc@10 while gcc@11 is already installed and
# linked, we want to overwrite those links and use gcc@10 instead.
brew link --overwrite $brew_package_name

# brew's llvm packages do not create symlinks as they might interfere with Apple's clang installation.
# Hence, we need to do it explicitly.
if [ "$package_name" == "llvm" ] && [[ $# -eq 2 ]]; then
    install_prefix=$(brew --prefix $brew_package_name)/bin
    ln -s $install_prefix/clang++ $install_prefix/clang++-$2 || true
    echo "$install_prefix" >> $GITHUB_PATH
fi
