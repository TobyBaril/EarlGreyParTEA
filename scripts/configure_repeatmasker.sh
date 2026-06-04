#!/usr/bin/env bash
set -euo pipefail

DB_PATH="$1"

LIB_DIR="${CONDA_PREFIX}/share/RepeatMasker/Libraries"
FAMDB_DIR="${LIB_DIR}/famdb/"

echo "Creating symbolic links from ${FAMDB_DIR} to Dfam DB at ${DB_PATH}"
for file in ${DB_PATH}/*.h5
do
    link=${FAMDB_DIR}/$(basename $file)
    ln -s $file $link
done

echo "Backing up mini init partition"
mv ${FAMDB_DIR}/min_init.0.h5 ${FAMDB_DIR}/min_init.0.h5.bak

echo "Running RepeatMasker configure..."
cd "${CONDA_PREFIX}/share/RepeatMasker"

perl ./configure \
    -libdir "$LIB_DIR" \
    -trf_prgm "${CONDA_PREFIX}/bin/trf" \
    -rmblast_dir "${CONDA_PREFIX}/bin" \
    -hmmer_dir "${CONDA_PREFIX}/bin" \
    -abblast_dir "${CONDA_PREFIX}/bin" \
    -crossmatch_dir "${CONDA_PREFIX}/bin" \
    -default_search_engine rmblast

echo "Marking configuration complete"
touch "$FAMDB_DIR/.earlgrey.config.complete"

echo "RepeatMasker configuration complete."
