#!/bin/bash

set -e

echo "Make sure cado submodule is initialized and updated"
cd ..
git submodule update --init

cd code/cado
git stash
CURRENT_COMMIT=$(git rev-parse --short HEAD)
CURRENT_COMMIT_DIR="build-${CURRENT_COMMIT}"
for patch_file in ../../patches/*; do
    patch -p1 < $patch_file
done
cd ..

MPI="0" CFLAGS="-O2 -DSUPPORT_LARGE_Q" CXXFLAGS="-O2 -DSUPPORT_LARGE_Q" make -f makefile.binaries
