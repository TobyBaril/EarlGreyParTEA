# Conda Release Automation

This directory contains the conda recipe for **earlgrey-partea**.

## How the release workflow operates

When a new GitHub Release is published (i.e. a tag in the format `vX.Y.Z` is released), the workflow defined in `.github/workflows/conda-release.yml` automatically:

1. **Extracts the version** from the release tag (stripping the leading `v`).
2. **Downloads the release tarball** and computes its `sha256` checksum.
3. **Updates `conda/meta.yaml`** in-place — replacing the version, `sha256`, and resetting the build number to `0`.
4. **Commits and pushes** the updated `meta.yaml` back to `main` with `[skip ci]` to prevent a CI loop.
5. **Lints** the conda recipe using `conda build --check`.
6. **Builds** the conda package using `conda build`.
7. **Uploads** the built package to the Anaconda channel [`toby_baril_bio`](https://anaconda.org/toby_baril_bio).

## Required repository secret

Before your first release you must add your Anaconda API token as a repository secret:

1. Go to **Settings → Secrets and Variables → Actions → New repository secret**
2. Name: `ANACONDA_API_TOKEN`
3. Value: your Anaconda.org API token (generate one at <https://anaconda.org/toby_baril_bio/settings/access>)

## Anaconda channel

Built packages are published to: <https://anaconda.org/toby_baril_bio>

Install the package with:

```bash
conda install -c toby_baril_bio -c conda-forge earlgrey-partea
```
