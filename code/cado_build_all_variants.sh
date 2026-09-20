#!/bin/bash

MPI_PATH="/usr/local/openmpi-5.0.8"

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

rm -rf $CURRENT_COMMIT_DIR
echo "Compiling all variants of cado for commit ${CURRENT_COMMIT_DIR}"
mkdir $CURRENT_COMMIT_DIR

if [ -d "build" ]; then
    echo "Backing up current build folder"
    rm -rf build.bkp
    mv build build.bkp
fi

build () {
    BUILD_NAME=$2
    FINAL_BUILD_DIR="${1}/${BUILD_NAME}"
    LOG_FILE="./${BUILD_NAME}"
    MPI="$3"
    FLAGS_SIZE="$4"
    echo "Compiling cado-nfs build '${BUILD_NAME}'"
    MPI="$MPI" FLAGS_SIZE="$FLAGS_SIZE" CFLAGS="-O2 -DSUPPORT_LARGE_Q" CXXFLAGS="-O2 -DSUPPORT_LARGE_Q" make -f makefile.binaries > "${LOG_FILE}.log" 2> "${LOG_FILE}.err"
    mv build "${FINAL_BUILD_DIR}"
    mv "${LOG_FILE}.log" "${LOG_FILE}.err" "${FINAL_BUILD_DIR}"
    echo "Finished building, stored in: ${FINAL_BUILD_DIR}"
}

build "$CURRENT_COMMIT_DIR" build-no-mpi-regular-lpb "0" ""
build "$CURRENT_COMMIT_DIR" build-with-mpi-regular-lpb "${MPI_PATH}" ""
build "$CURRENT_COMMIT_DIR" build-no-mpi-large-lpb "0" "-DSIZEOF_P_R_VALUES=8 -DSIZEOF_INDEX=8"
build "$CURRENT_COMMIT_DIR" build-with-mpi-large-lpb "${MPI_PATH}" "-DSIZEOF_P_R_VALUES=8 -DSIZEOF_INDEX=8"

if [ -d "build.bkp" ]; then
    echo "Restoring previous build folder"
    mv build.bkp build
fi

