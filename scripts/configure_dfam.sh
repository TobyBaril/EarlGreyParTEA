#!/usr/bin/env bash
# configure_dfam.sh — Dfam library setup guide for EarlGrey ParTEA
#
# EarlGrey 7.3.0+ (Dfam 4.0 / FamDB 3.0.0)
# ─────────────────────────────────────────
# RepeatMasker is pre-configured by the conda post-install hook.
# No 'perl ./configure' step is required.
#
# To download the Dfam 4.0 library partitions, run the interactive
# download tool provided by the famdb conda package:
#
#   download_dfam.py
#
# This will guide you through selecting which partitions to download
# (curated consensus, HMMs, uncurated sequences, etc.) and writes the
# HDF5 files to the correct famdb Libraries directory automatically.
#
# ─────────────────────────────────────────
# Legacy note (EarlGrey ≤7.2, Dfam 3.9)
# ─────────────────────────────────────────
# If you are using an older EarlGrey environment (< 7.3.0) with Dfam 3.9,
# the original manual setup steps were:
#
#   cd $CONDA_PREFIX/share/RepeatMasker/Libraries/famdb/
#   curl -o 'dfam39_full.#1.h5.gz' \
#     'https://dfam.org/releases/current/families/FamDB/dfam39_full.[0-16].h5.gz'
#   gunzip -f *.gz
#   mv min_init.0.h5 min_init.0.h5.bak
#   cd $CONDA_PREFIX/share/RepeatMasker/
#   perl ./configure \
#       -libdir $CONDA_PREFIX/share/RepeatMasker/Libraries \
#       -trf_prgm $CONDA_PREFIX/bin/trf \
#       -rmblast_dir $CONDA_PREFIX/bin \
#       -hmmer_dir $CONDA_PREFIX/bin \
#       -default_search_engine rmblast
#   touch $CONDA_PREFIX/share/RepeatMasker/Libraries/famdb/.earlgrey.config.complete
