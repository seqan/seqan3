#!/usr/bin/env bash
set -x
set -e

echo 'APT::Acquire::Retries "5";' | sudo tee -a /etc/apt/apt.conf.d/80-retries > /dev/null
sudo add-apt-repository --no-update --yes ppa:ubuntu-toolchain-r/ppa
sudo add-apt-repository --no-update --yes ppa:ubuntu-toolchain-r/test
sudo apt-get update
