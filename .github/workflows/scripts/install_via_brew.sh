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

brew list --versions $package_name && brew unlink $package_name || true
brew list --versions $package_name && brew upgrade --force-bottle $brew_package_name || brew install --force-bottle $brew_package_name
brew link --overwrite $brew_package_name

