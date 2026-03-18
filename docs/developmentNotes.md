# Development record following initial release v0.1.2

### Implementation of new feature: Number of TE families added as genome number increases

**Goal:** Visualise how the pangenome TE library saturates as additional genomes are
included in the pipeline — a saturation/accumulation curve analogous to pangenome
open/closed curves used in bacterial genomics.

**Approach:** Parse the `.clstr` membership file produced automatically by `cd-hit-est`
during the `cluster_all_species` step (`combined_all_species.clstrd.fa.clstr`). No
re-clustering is performed. Each cluster in the file represents one unique TE family
in the final library, and every member sequence carries a species-prefix header
(`{species}_`, `REPMASKER_{repspec}_`, or `CUSTOM_`) that records its origin.

**Algorithm:**
1. Parse `.clstr` → map of cluster ID to the set of source species/libraries that
   contributed at least one sequence to that cluster.
2. Assign clusters containing any `REPMASKER_` or `CUSTOM_` sequence as the x=0
   baseline, representing families already known before any new genome is added.
   Because clustering is run on all sequences together, any de novo family that is
   redundant with an existing library family is absorbed into the same cluster and
   is therefore NOT counted as a novel addition — the correct behaviour.
3. For N random permutations of genome addition order, compute the cumulative number
   of unique clusters (TE families) discovered at each step — a cluster is first
   attributed to the earliest genome in the permutation order that contributed a
   member sequence.
4. Report mean ± 95% percentile CI across all permutations as the saturation curve.

**Existing library handling:**
- If `repeatmasker_species` is set, RepeatMasker library sequences are prefixed
  `REPMASKER_{repspec}_` during clustering and contribute to the x=0 baseline.
- If `custom_library` is set, custom library sequences are prefixed `CUSTOM_` and
  similarly contribute to the x=0 baseline.
- If neither is set, the x axis starts at 1 (first genome) and the y intercept is
  omitted since it is trivially zero.
- If `skip_clustering: True`, no `.clstr` file is produced by cd-hit-est. An empty
  sentinel `.clstr` is written so Snakemake dependency tracking still works. The
  saturation script detects this (file size == 0) and falls back to counting raw
  sequences per per-genome strained FASTA file with a warning. In this mode,
  cross-genome redundancy cannot be accounted for.

**Outputs written to `{outdir}/combinedLibraries/`:**
- `saturation_plot.pdf` — saturation curve (mean line + shaded 95% CI band)
- `saturation_data.tsv` — columns: `n_genomes, mean_unique_families, ci_lower_95,
  ci_upper_95`

**New files added to the pipeline:**
- `scripts/saturation_plot.py` — parsing, permutation sampling, and plotting logic
- `rules/saturation.smk` — Snakemake rule triggered after `cluster_all_species`

**Changes to existing files:**
- `rules/clustering.smk` — expose `combined_all_species.clstrd.fa.clstr` as a
  declared output of `cluster_all_species`; touch an empty sentinel when
  `skip_clustering: True`
- `Snakefile` — include `saturation.smk` for "full" and "libconstruct" modes;
  add saturation outputs to `rule all` targets for those modes
- `config/config.yaml` — add `saturation_permutations: 100` (tunable)

**Design rationale:**
- Using the final clustered `.clstr` as a proxy avoids the O(N) cost of re-running
  cd-hit-est on every genome subset. This is the standard approach in pangenome
  accumulation analyses (Tettelin et al. 2005 framework) and gives accurate
  attribution because all sources were clustered together in a single run.
- The permutation approach (rather than exhaustive enumeration of all N! orderings)
  is standard; 100 permutations converges well for datasets of up to ~50 genomes.
- Species name prefixes in headers are matched longest-first to handle names that
  are substrings of other names.

**Verification checklist:**
- [ ] Run pipeline on 3+ genomes; confirm `.clstr` file exists after clustering step
- [ ] Run saturation script manually; confirm TSV has monotonically non-decreasing
      mean and the correct number of row steps
- [ ] With `repeatmasker_species` set: confirm x=0 baseline > 0 in plot
- [ ] With `custom_library` set: confirm same baseline behaviour
- [ ] With `skip_clustering: True`: confirm graceful fallback and warning in log
- [ ] With a single genome: curve is a single point (CI width = 0)

**To Test:**
Make conda environment from development branch:
```bash
cd /data/toby/EarlGreyParTEA
conda-build purge-all
conda build conda/
conda create -n test_saturation --use-local earlgrey-partea
conda activate test_saturation
```

symlink earlgrey libraries and configure RepeatMasker:
```bash
ln -s /data/toby/tools/earlgrey_databases/Libraries/famdb/* /data/toby/miniforge3/envs/test_saturation/share/RepeatMasker/Libraries/famdb/
ln -s /data/toby/tools/earlgrey_databases/Libraries/RMRB.embl /data/toby/miniforge3/envs/test_saturation/share/RepeatMasker/Libraries/RMRB.embl
ln -s /data/toby/tools/earlgrey_databases/Libraries/RMRBSeqs.embl /data/toby/miniforge3/envs/test_saturation/share/RepeatMasker/Libraries/RMRBSeqs.embl

cd /data/toby/miniforge3/envs/test_saturation/share/RepeatMasker/
/data/toby/miniforge3/envs/test_saturation/bin/perl ./configure
cd /data/toby/EarlGreyParTEA
```

Make a test directory with 4 Z. tritici genomes and a config only doing de novo library construction:
```bash
mkdir -p /data/toby/testDIR/saturationTests/1_4genomesNoRepMasker
cd /data/toby/testDIR/saturationTests/1_4genomesNoRepMasker
# copy some Z. tritici genomes
cp /legserv/NGS_data/Zymoseptoria/Zt_Reference_genomes/19Pangenome_genomes/IPO323/Zymoseptoria_tritici.MG2.dna.toplevel.mt+.fa /data/toby/testDIR/saturationTests/1_4genomesNoRepMasker/IPO323.fa

cp /legserv/NGS_data/Zymoseptoria/Zt_Reference_genomes/19Pangenome_genomes/1A5/ST99CH_1A5.fa /data/toby/testDIR/saturationTests/1_4genomesNoRepMasker/1A5.fa

cp /legserv/NGS_data/Zymoseptoria/Zt_Reference_genomes/19Pangenome_genomes/YEQ92/YEQ92.fa /data/toby/testDIR/saturationTests/1_4genomesNoRepMasker/YEQ92.fa

cp /legserv/NGS_data/Zymoseptoria/Zt_Reference_genomes/19Pangenome_genomes/Aus01/Aus01.fa /data/toby/testDIR/saturationTests/1_4genomesNoRepMasker/Aus01.fa
```

Create a config.yaml using the earlGreyParTEA command:
```bash
earlGreyParTEA_LibConstruct --generate-config 1_4genomesNoRepMasker_config.yaml
```

Run the pipeline in libconstruct mode:
```bash
earlGreyParTEA_LibConstruct \
    -c /data/toby/testDIR/saturationTests/1_4genomesNoRepMasker/1_4genomesNoRepMasker_config.yaml \
    -t 8
```

We have a saturation plot and data! Check that the curve is monotonically non-decreasing and that the x=0 baseline is zero (since no RepeatMasker library was used). Then repeat with `repeatmasker_species: "Zymoseptoria"` set in the config to confirm a non-zero baseline. Finally, test with `skip_clustering: True` to confirm the fallback logic works and a warning is issued in the log.

Create a config.yaml using the earlGreyParTEA command:
```bash
cd /data/toby/testDIR/saturationTests
earlGreyParTEA_LibConstruct --generate-config 2_4genomesWithRepMasker_config.yaml
```

Run the pipeline in libconstruct mode:
```bash
earlGreyParTEA_LibConstruct \
    -c /data/toby/testDIR/saturationTests/2_4genomesWithRepMasker_config.yaml \
    -t 8
```

### Issues with config
Daniel had an issue where he used a custom config, and it was still pulling some information from the `config/config.yaml` file. This should not be used at all. The config file generated by the command should be fully self-contained and not rely on any values from the default config. This is to ensure reproducibility and that users can share their config files without worrying about hidden dependencies on the default config.

**Root cause:** The `Snakefile` contained the directive `configfile: "config/config.yaml"`. Snakemake processes `configfile:` directives before merging any CLI `--configfile` argument, so `config/config.yaml` was always loaded first. Any key present in the default config that was absent from the user's config would silently bleed through into the merged `config` dict seen by the workflow — a hidden dependency invisible to the user.

**Fix:**
- `Snakefile` — removed the `configfile: "config/config.yaml"` directive entirely. The workflow now only uses the config explicitly passed by the user via `--configfile`. No values from the default config can bleed through.
- `scripts/on_start_functions.py` — added `saturation_permutations` (default `100`) to the centralised `defaults` dict in `validate_parameters`. This was the one parameter that had previously relied on scattered `config.get(..., 100)` fallbacks rather than the single authoritative defaults dict that all other optional parameters use.

All other optional parameters already had Python-level defaults in `validate_parameters`, so removing the `configfile:` directive is safe. The `--generate-config` command already produces fully self-contained configs with every field explicitly set.

### RepeatMasker general library cache race condition

When multiple RepeatMasker jobs start simultaneously in a freshly configured environment, they all attempt to build the `Libraries/general` cache at the same time. The first process to finish writes the cache; every other process encounters a partially-written `general.working` directory and fails with:

```
RepeatMasker::createLib(): Error invoking .../makeblastdb on file .../general.working/is.lib.
```

This only happens on the very first RepeatMasker run in a new conda environment. Subsequent runs succeed because the cache already exists. Re-submitting the pipeline after the failure resolves it, but that requires manual intervention.

**Fix:** A new `repeatmasker_warmup` rule in `rules/lib_construct.smk` runs a single RepeatMasker job on an inline dummy FASTA before any genome jobs start. It checks whether the `Libraries/general` cache directory already exists (conda path: `$(which RepeatMasker | sed 's|/bin/RepeatMasker$|/share/RepeatMasker/Libraries/general|')`); if it does, it exits immediately. If not, it runs RepeatMasker with `-lib dummy.fa` on the dummy sequence to trigger the cache build. The output is a sentinel file (`{outdir}/.repeatmasker_cache_ready`, created via Snakemake's `touch()`). Every RepeatMasker rule lists this sentinel as an input (`_cache`), ensuring the warmup completes before any parallel genome jobs are dispatched. On subsequent runs the sentinel exists and Snakemake skips the warmup entirely.

`rules/annotate_simple.smk` was also removed at this point — it was a dead file not referenced by the `Snakefile` or any other rule file.

**Files changed:**
- `rules/lib_construct.smk` — added `repeatmasker_warmup` rule; added `_cache` sentinel input to `repeatmasker` and `repeatmasker_custom` rules
- `rules/annotate.smk` — added `_cache` sentinel input to `repeatmasker_annotation` rule
- `rules/annotate_simple.smk` — deleted (dead file)

## Test Release v0.1.3
- [x] Bump version to 0.1.3 in all relevant files:
  - `earlGreyParTEA`
  - `earlGreyParTEA_AnnotationOnly`
  - `earlGreyParTEA_LibConstruct`
  - `conda/meta.yaml`
- [x] Build conda package and test install in a fresh environment
- [x] Run pipeline on test dataset with new config parameters; confirm expected outputs and behaviour

Build conda package:
```bash
cd /data/toby/EarlGreyParTEA
conda build conda/
conda create -n test_013 --use-local earlgrey-partea
conda activate test_013
```

Configure RepeatMasker in the new environment (if not already done):
```bash
ln -s /data/toby/tools/earlgrey_databases/Libraries/famdb/* /data/toby/miniforge3/envs/test_013/share/RepeatMasker/Libraries/famdb/
ln -s /data/toby/tools/earlgrey_databases/Libraries/RMRB.embl /data/toby/miniforge3/envs/test_013/share/RepeatMasker/Libraries/RMRB.embl
ln -s /data/toby/tools/earlgrey_databases/Libraries/RMRBSeqs.embl /data/toby/miniforge3/envs/test_013/share/RepeatMasker/Libraries/RMRBSeqs.embl
cd /data/toby/miniforge3/envs/test_013/share/RepeatMasker/
/data/toby/miniforge3/envs/test_013/bin/perl ./configure
cd /data/toby/EarlGreyParTEA
```

Run a full test with the new config parameters:
```bash
cd /data/toby/testDIR/
# make a config with the new parameters set
earlGreyParTEA --generate-config 3_updated0.1.3.config.yaml

# run the full pipeline with the new config
earlGreyParTEA \
    -c /data/toby/testDIR/3_updated0.1.3.config.yaml \
    -t 32
```

deactivate the environment, delete it, purge the build, then build again to test the RepeatMasker cache issues
```bash
conda deactivate
conda remove -n test_013 --all
conda-build purge-all

cd /data/toby/EarlGreyParTEA
conda build conda/
conda create -n test_013 --use-local earlgrey-partea
conda activate test_013

ln -s /data/toby/tools/earlgrey_databases/Libraries/famdb/* /data/toby/miniforge3/envs/test_013/share/RepeatMasker/Libraries/famdb/
ln -s /data/toby/tools/earlgrey_databases/Libraries/RMRB.embl /data/toby/miniforge3/envs/test_013/share/RepeatMasker/Libraries/RMRB.embl
ln -s /data/toby/tools/earlgrey_databases/Libraries/RMRBSeqs.embl /data/toby/miniforge3/envs/test_013/share/RepeatMasker/Libraries/RMRBSeqs.embl
cd /data/toby/miniforge3/envs/test_013/share/RepeatMasker/
/data/toby/miniforge3/envs/test_013/bin/perl ./configure
cd /data/toby/EarlGreyParTEA

# run the full pipeline with the new config
earlGreyParTEA \
    -c /data/toby/testDIR/3_updated0.1.3.config.yaml \
    -t 32
```

This works without any manual intervention, confirming the warmup logic correctly handles the RepeatMasker cache race condition. Check the logs to confirm the warmup rule ran on the first run and was skipped on the second run. Also check that the saturation plot and data were generated correctly in both runs.

I will commit these changes to the `development` branch and then merge into `main` for release v0.1.3.

### Relative path failure when pipeline is run from a different working directory

Rules in `lib_construct.smk` and `annotate.smk` that `cd` into a subdirectory (e.g. `build_db`, `repeatmodeler`, `testrainer`, `heliano_detection`) used input/param values such as `{input.masked}` or `{input.genome}` **after** the `cd`. When these paths were relative (e.g. `output_dir: condaPull`), the shell could no longer resolve them from the new working directory, producing errors such as:

```
Command line fasta file condaPull/genome2_EarlGrey/genome2.prep does not exist!
FileNotFoundError: No such file or directory: '/data/toby/testDIR/condaPull/genome2_EarlGrey/genome2_heliano/condaPull/genome2_EarlGrey/genome2.prep'
```

**Root cause:** Snakemake expands wildcard patterns using whatever string the user put in `output_dir`, `genome`, etc. If those strings are relative, all derived paths are relative. Any shell block that runs `cd` before referencing an input path then silently looks in the wrong place.

**Fix:** In `validate_parameters` (`scripts/on_start_functions.py`), immediately before the output-directory setup block, all user-supplied file/directory paths are converted to absolute paths with `os.path.abspath`:

```python
config['output_dir'] = os.path.abspath(config['output_dir'])
config['genome'] = {sp: os.path.abspath(p) for sp, p in config['genome'].items()}
if config.get('custom_library'):
    config['custom_library'] = os.path.abspath(config['custom_library'])
if config.get('annotation_library'):
    config['annotation_library'] = os.path.abspath(config['annotation_library'])
```

Because `validate_parameters` is called at parse time (before any rules run), all wildcard expansions of `{outdir}` and all input genome paths are already absolute by the time Snakemake builds the DAG. No per-rule `$(realpath ...)` patches are needed.

**Files changed:**
- `scripts/on_start_functions.py` — added path absolutization block in `validate_parameters`

## Build and upload to toby_baril_bio channel on Anaconda Cloud

```bash
cd /data/toby/EarlGreyParTEA
conda build conda/
anaconda login
anaconda upload /data/toby/miniforge3/conda-bld/noarch/earlgrey-partea-0.1.3-py_0.conda
```

## Release 0.1.5 Feature Updates

### Extended config generation: `--genome-dir` and `--from-csv`

**Goal:** Lower the barrier for new users by letting them auto-populate the
`genome:` and `species:` blocks of a config file from either a directory of
FASTA files or a simple CSV spreadsheet, rather than editing the template by
hand.

**Three supported workflows:**
1. **Blank template** (unchanged) — `--generate-config my_config.yaml` produces
   a fully populated config with placeholder genome paths. The user manually
   fills in paths.
2. **From a genome directory** — `--generate-config my_config.yaml --genome-dir
   /path/to/genomes/ [--output-dir /path/to/results/]` scans the directory for
   FASTA files (`.fa`, `.fasta`, `.fna`, `.fa.gz`, `.fasta.gz`), uses the
   filename-without-extension as the species name, and auto-populates the
   `genome:` and `species:` blocks with absolute paths.
3. **From a CSV file** — `--generate-config my_config.yaml --from-csv
   genomes.csv [--output-dir /path/to/results/]` reads a CSV with at minimum
   two columns (`species` and `genome_path`) and populates the config
   accordingly. Relative genome paths in the CSV are resolved relative to the
   CSV file's own location.

**Important design decision:** `output_dir` is always a single value supplied
via `--output-dir` on the command line. It is never read from the CSV — keeping
a clean separation between the per-genome input metadata (CSV) and the
run-level output location.

**Species name sanitization:** Species names are passed through a sanitizer
that replaces characters unsafe in YAML keys or TE sequence headers with `_`.
Names beginning with a digit are prefixed `s_`. Both transformations emit a
warning so users are aware of the change.

**Implementation:**
- `scripts/generate_config.py` — new Python helper containing all genome-dir
  scanning and CSV parsing logic. Accepts `--output`, `--mode`
  (`full|libconstruct|annotate`), `--genome-dir`, `--from-csv`, and
  `--output-dir`. Each bash entry point passes its own mode string, keeping
  the generated `pipeline_mode:` field correct for the command used.
- `earlGreyParTEA`, `earlGreyParTEA_LibConstruct`,
  `earlGreyParTEA_AnnotationOnly` — new optional flags `--genome-dir`,
  `--from-csv`, `--output-dir` added to argument parsing. When either
  `--genome-dir` or `--from-csv` is present, the script locates
  `generate_config.py` using the same search-path logic used for the
  `Snakefile`, then delegates to Python. Without either flag, the existing
  inline blank-template behaviour is unchanged (backward compatible, no Python
  dependency for the plain template path).
- `--genome-dir` and `--from-csv` are mutually exclusive; supplying both
  produces an error.

**Files changed:**
- `scripts/generate_config.py` — new file
- `earlGreyParTEA` — updated usage, argument parsing, and config-generation
  dispatch block
- `earlGreyParTEA_LibConstruct` — same
- `earlGreyParTEA_AnnotationOnly` — same

