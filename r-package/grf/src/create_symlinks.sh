#!/usr/bin/env bash
# Symlinks the grf C++ dependencies.
#
# This only needs to be run once, from within this src directory:
#    `./create_symlinks.sh`

ln -s ../../../core/src/
ln -s ../../../core/third_party/optional/
ln -s ../bindings/* .
