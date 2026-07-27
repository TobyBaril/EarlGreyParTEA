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

### SLURM cluster submission via `--slurm`

**Goal:** When `--slurm` is passed, submit each Snakemake rule as an individual
`sbatch` job to a SLURM cluster, with CPUs, memory, and runtime derived from the
`resources:` block of each rule. The user specifies the partition and optionally
account/extra flags; sensible conservative defaults are baked into each rule so
jobs are likely to complete without manual tuning.

---

#### Approach: Snakemake 8+ native SLURM executor plugin

The environment runs Snakemake 9.6.2. The correct approach is to use the official
`snakemake-executor-plugin-slurm` package, which is installed alongside Snakemake
and activated with `--executor slurm`. This is cleaner than the legacy
`--cluster "sbatch ..."` mode: job tracking, cancellation, and status polling
are handled internally by the plugin — no custom `slurm_status.py` script is needed.

Per-rule resource control flows through the existing `threads:` and `resources:`
directives. The plugin maps `mem_mb` → `--mem`, `runtime` → `--time` (minutes),
`threads` → `--cpus-per-task`, `slurm_partition` → `--partition`, and
`slurm_account` → `--account`.

`snakemake-executor-plugin-slurm` must be added as a run dependency in
`conda/meta.yaml`.

---

#### New CLI flags (all three entry points)

| Flag | Description |
|------|-------------|
| `--slurm` | Enable SLURM submission mode |
| `--slurm-jobs N` | Max jobs submitted concurrently (default: total number of genomes × 3, capped at 200) |
| `--slurm-partition PART` | SLURM partition/queue to submit to (required with `--slurm`) |
| `--slurm-account ACCT` | SLURM account string (optional, omitted from `sbatch` if empty) |
| `--slurm-extra "FLAGS"` | Any additional raw `sbatch` flags (e.g. `"--constraint=avx2"`) |

`-t/--threads` changes meaning in SLURM mode: it becomes **CPUs per job** rather
than total cores across the machine. Concurrency is controlled by `--slurm-jobs`
instead.

---

#### Snakemake invocation in SLURM mode

```bash
# Constructed in the bash wrapper:
DEFAULT_RESOURCES=("slurm_partition=$SLURM_PARTITION")
[ -n "$SLURM_ACCOUNT" ] && DEFAULT_RESOURCES+=("slurm_account=$SLURM_ACCOUNT")
[ -n "$SLURM_EXTRA" ]   && DEFAULT_RESOURCES+=("slurm_extra=$SLURM_EXTRA")

snakemake \
    --snakefile "$SNAKEFILE" \
    --configfile "$CONFIG" \
    --config slurm_mode=true \
    --executor slurm \
    --default-resources "${DEFAULT_RESOURCES[@]}" \
    --jobs "$SLURM_JOBS" \
    --cores "$THREADS" \
    --latency-wait 60 \
    --retries 1 \
    $DRY_RUN $UNLOCK $RERUN
```

Key flags:
- `--executor slurm` — activates the `snakemake-executor-plugin-slurm` plugin.
  Job submission, tracking, and cancellation are handled automatically.
- `--default-resources` — sets partition (and optionally account/extra flags)
  for all rules that do not specify them in their own `resources:` block.
- `--config slurm_mode=true` — injects a flag into the Snakemake config that the
  thread lambdas read to skip the `// len(SPECIES_LIST)` division (see below).
- `--latency-wait 60` — gives network filesystems 60 s to make output files visible
  after a job completes. Critical on NFS-mounted HPC scratch.
- `--retries 1` — one automatic retry at 2× memory on failure (via `attempt` scaling).
- No `--cores` is passed in SLURM mode; concurrency is governed by `--jobs` alone.

---

#### Thread lambda modification

Current lambdas divide `workflow.cores` by the number of species so all per-genome
jobs share the machine. In SLURM mode each job has its own allocation, so the
division must not happen. The lambdas in `lib_construct.smk` and `annotate.smk`
are updated to check `config.get("slurm_mode", False)`:

```python
# Before (local only)
threads: lambda wildcards: max(1, min(workflow.cores // len(SPECIES_LIST), 64))

# After (mode-aware)
threads: lambda wildcards: (
    max(1, min(workflow.cores, 64))
    if config.get("slurm_mode", False)
    else max(1, min(workflow.cores // len(SPECIES_LIST), 64))
)
```

The `cluster_all_species` and saturation rules that already use `workflow.cores`
directly (not divided) do not need modification.

---

#### Resource requirements per rule

`runtime` is specified in **minutes** (SLURM accepts a bare integer as minutes).
`mem_mb` scales with retry attempt (`attempt`) so transient OOM failures trigger an
automatic retry at higher memory before failing permanently. Rules that currently
lack a `resources:` block gain one.

| Rule | threads (SLURM) | mem_mb | runtime (min) | Notes |
|------|-----------------|--------|---------------|-------|
| `repeatmasker_warmup` | 1 | 4 000 | 30 | Cache build; very short |
| `prep_genome` | 1 | 4 000 | 60 | Fast I/O-bound step |
| `extract_repeatmasker_library` | 1 | 4 000 | 30 | famdb query; fast |
| `build_db` | 1 | 8 000 × attempt | 120 | BuildDatabase |
| `repeatmasker` | `--threads` value | 16 000 × attempt | 10 080 (1 week) | Initial masking |
| `repeatmasker_custom` | `--threads` value | 16 000 × attempt | 10 080 (1 week) | Custom lib masking |
| `repeatmodeler` | `--threads` value | 32 000 × attempt | 10 080 (1 week) | Most expensive step |
| `testrainer` | `--threads` value | 16 000 × attempt | 10 080 (1 week) | BLAST cycles |
| `cluster_all_species` | up to 32 | 32 000 × attempt | 480 (8 h) | cd-hit; **locally executed** (see below) |
| `repeatmasker_annotation` | `--threads` value | 16 000 × attempt | 480 (8 h) | |
| `heliano_detection` | `--threads` value | 8 000 × attempt | 480 (8 h) | |
| `merge_repeats` | up to 16 | 8 000 × attempt | 240 (4 h) | |
| `generate_summary_charts` | 1 | 4 000 | 30 | R/shell; fast |
| `calculate_divergence` | up to 16 | 8 000 × attempt | 480 (8 h) | |
| `sweep_up_files` | 1 | 2 000 | 15 | File copies only |
| `generate_softmasked_genome` | 1 | 8 000 | 60 | bedtools |
| `saturation_plot` | 1 | 4 000 | 30 | Python; **locally executed** (see below) |

Default retry count (`--retries 1` passed to Snakemake in SLURM mode) means a
rule gets one automatic retry at 2× memory before Snakemake marks it failed.

---

#### Locally-executed rules in cluster mode

In Snakemake 8+, the SLURM executor plugin re-invokes Snakemake on the cluster
node to execute each rule, including rules with `run:` Python blocks and `script:`
directives. This means **all rules — including `cluster_all_species` and
`saturation_plot` — are submitted as proper SLURM jobs** and no longer run on the
submit node. The memory and runtime resources set for those rules are respected.

This resolves the Snakemake 7 limitation and means no special handling or user
warnings are needed for those rules.

---

#### Config additions

Three new optional keys added to `config.yaml` and the `generate_config.py`
templates. They are only used when `--slurm` is passed:

```yaml
# SLURM cluster settings (only used with --slurm flag)
slurm_partition: "long"   # partition/queue to submit to
slurm_account: ""         # account string (leave empty if not required)
slurm_extra: ""           # any extra sbatch flags, e.g. "--constraint=avx2"
```

These are documented in the config but do not affect local runs.

---

#### Files to create/modify

- `conda/meta.yaml` — add `snakemake-executor-plugin-slurm` as a run dependency
- `earlGreyParTEA` — add `--slurm`, `--slurm-jobs`, `--slurm-partition`,
  `--slurm-account`, `--slurm-extra` flags; build `--default-resources` string;
  invoke Snakemake with `--executor slurm`
- `earlGreyParTEA_LibConstruct` — same
- `earlGreyParTEA_AnnotationOnly` — same
- `rules/lib_construct.smk` — update thread lambdas for SLURM mode; add `runtime`
  and `mem_mb` resources to all rules currently missing them
- `rules/annotate.smk` — same
- `rules/clustering.smk` — add `runtime` resource to `cluster_all_species`
- `rules/saturation.smk` — add `runtime` and `mem_mb` resources
- `config/config.yaml` — add SLURM config keys
- `scripts/generate_config.py` — add SLURM keys to all three template strings

---

#### Verification checklist

- [ ] Dry run with `--slurm` shows correct `sbatch` command printed per job
- [ ] Submit one genome in SLURM mode; confirm job appears in `squeue` with correct
      CPUs, memory, and time limit
- [ ] A deliberately OOM job retries at 2× memory automatically (attempt-scaled mem_mb)
- [ ] `cluster_all_species` is submitted as a proper SLURM job (not run on submit node)
- [ ] Cancelling Snakemake with Ctrl-C cancels all submitted jobs via the plugin
- [ ] `--latency-wait 60` prevents false "output file missing" failures on NFS scratch
- [ ] Local mode behaviour is completely unchanged (no regression)

## Test Release v0.1.5

### Version bump checklist

Before building, update the version string in all relevant files:

- [ ] `earlGreyParTEA` — `VERSION="0.1.5"`
- [ ] `earlGreyParTEA_LibConstruct` — `VERSION="0.1.5"`
- [ ] `earlGreyParTEA_AnnotationOnly` — `VERSION="0.1.5"`
- [ ] `conda/meta.yaml` — `version: "0.1.5"` and update `sha256` after tagging

### Build and install

```bash
cd /data/toby/EarlGreyParTEA
conda-build purge-all
conda build conda/
conda create -n test_015 --use-local earlgrey-partea
conda activate test_015
```

Configure RepeatMasker in the new environment:

```bash
ln -s /data/toby/tools/earlgrey_databases/Libraries/famdb/* /data/toby/miniforge3/envs/test_015/share/RepeatMasker/Libraries/famdb/
ln -s /data/toby/tools/earlgrey_databases/Libraries/RMRB.embl /data/toby/miniforge3/envs/test_015/share/RepeatMasker/Libraries/RMRB.embl
ln -s /data/toby/tools/earlgrey_databases/Libraries/RMRBSeqs.embl /data/toby/miniforge3/envs/test_015/share/RepeatMasker/Libraries/RMRBSeqs.embl
cd /data/toby/miniforge3/envs/test_015/share/RepeatMasker/
/data/toby/miniforge3/envs/test_015/bin/perl ./configure
cd /data/toby/EarlGreyParTEA
```

---

### Test 1: `--genome-dir` config generation

Scan the existing test genome directory and auto-populate a config. The four
Z. tritici FASTAs in `1_4genomesNoRepMasker/` serve as a ready-made test set.

```bash
mkdir -p /data/toby/testDIR/saturationTests/3_genomeDirTest
cd /data/toby/testDIR/saturationTests/3_genomeDirTest

earlGreyParTEA_LibConstruct \
    --generate-config 3_genomedir_config.yaml \
    --genome-dir /data/toby/testDIR/saturationTests/1_4genomesNoRepMasker/ \
    --output-dir /data/toby/testDIR/saturationTests/3_genomeDirTest/results

cat 3_genomedir_config.yaml
```

Confirm:
- All four genomes (IPO323, 1A5, YEQ92, Aus01) appear under `genome:` with absolute paths
- `species:` list contains all four names
- `output_dir:` is set to the value passed via `--output-dir`
- `pipeline_mode: "libconstruct"` (not `full`)
- The species name `1A5` is sanitized to `s_1A5` and a warning was printed

---

### Test 2: `--from-csv` config generation

Test `species` and `genome_path` columns. Include a species name starting with a digit (1A5) to
exercise the sanitizer, and one with a special character (e.g. a dot) to confirm replacement with `_`.

```bash
mkdir -p /data/toby/testDIR/saturationTests/4_fromCSVTest
cd /data/toby/testDIR/saturationTests/4_fromCSVTest

cat > genomes.csv << 'EOF'
species,genome_path
IPO323,/data/toby/testDIR/saturationTests/1_4genomesNoRepMasker/IPO323.fa
1A5,/data/toby/testDIR/saturationTests/1_4genomesNoRepMasker/1A5.fa
YEQ92,/data/toby/testDIR/saturationTests/1_4genomesNoRepMasker/YEQ92.fa
Aus01,/data/toby/testDIR/saturationTests/1_4genomesNoRepMasker/Aus01.fa
EOF

earlGreyParTEA_LibConstruct \
    --generate-config 4_fromcsv_config.yaml \
    --from-csv genomes.csv \
    --output-dir /data/toby/testDIR/saturationTests/4_fromCSVTest/results

cat 4_fromcsv_config.yaml
```

Confirm:
- `genome:` and `species:` are identical to Test 1 above
- `1A5` → `s_1A5` with a printed warning
- `output_dir:` matches `--output-dir`
- Specifying both `--genome-dir` and `--from-csv` at the same time produces an error:

```bash
earlGreyParTEA_LibConstruct \
    --generate-config should_fail.yaml \
    --genome-dir /data/toby/testDIR/saturationTests/1_4genomesNoRepMasker/ \
    --from-csv genomes.csv
# Expected: [ERROR] --genome-dir and --from-csv are mutually exclusive
```

---

### Test 3: SLURM dry run

A dry run should print the Snakemake execution plan and confirm that `--executor slurm` is
present in the constructed command, without submitting any real jobs.

```bash
cd /data/toby/testDIR/saturationTests

earlGreyParTEA_LibConstruct \
    -c /data/toby/testDIR/saturationTests/1_4genomesNoRepMasker/1_4genomesNoRepMasker_config.yaml \
    -t 8 \
    --slurm \
    --slurm-partition normal.1000h \
    --dry-run
```

Confirm:
- The printed `Executing:` line contains `--executor slurm`
- The printed line contains `--latency-wait 60` and `--retries 1`
- `--jobs` is set to `12` (4 genomes × 3, < 200 cap)
- Rules that were already complete are skipped; new rules show the expected job graph

Test that omitting `--slurm-partition` with `--slurm` fails immediately:

```bash
earlGreyParTEA_LibConstruct \
    -c /data/toby/testDIR/saturationTests/1_4genomesNoRepMasker/1_4genomesNoRepMasker_config.yaml \
    -t 8 \
    --slurm \
    --dry-run
# Expected: [ERROR] --slurm-partition is required when using --slurm
```

---

### Test 4: SLURM live submission

Replace `<your-partition>` with a real partition available on the cluster (check with `sinfo`).

```bash
cd /data/toby/testDIR/saturationTests

mkdir -p /data/toby/testDIR/saturationTests/5_slurmTest

earlGreyParTEA_LibConstruct \
    -c /data/toby/testDIR/saturationTests/1_4genomesNoRepMasker/1_4genomesNoRepMasker_config.yaml \
    -t 32 \
    --slurm \
    --slurm-partition normal.1000h \
    --slurm-jobs 12

# While running, check that jobs appear with correct CPUs, memory, and time:
squeue -u $USER -o "%.18i %.9P %.30j %.8u %.8T %.10M %.9l %.6C %.10m"
```

Confirm:
- Jobs appear in `squeue` with the correct partition
- CPUs match `-t 8`
- Memory and time limits match the per-rule values in the resource table (e.g. `repeatmodeler` should
  show `~32GB` and `7-00:00:00`)
- `cluster_all_species` and `saturation_plot` are submitted as proper SLURM jobs (not run on the
  submit node)
- Upon pipeline completion, `saturation_plot.pdf` and `saturation_data.tsv` exist in
  `combinedLibraries/`

Test `--slurm-account` flag:

```bash
earlGreyParTEA_LibConstruct \
    -c /data/toby/testDIR/saturationTests/1_4genomesNoRepMasker/1_4genomesNoRepMasker_config.yaml \
    -t 8 \
    --slurm \
    --slurm-partition normal.1000h \
    --slurm-account toby \
    --dry-run
# Confirm --default-resources contains slurm_account=<your-account> in the printed command
```

---

### Test 5: Local mode regression

Confirm that the SLURM flags have not broken standard local execution:

```bash
# Dry run in local mode — confirm no slurm flags appear in the Executing: line
earlGreyParTEA_LibConstruct \
    -c /data/toby/testDIR/saturationTests/1_4genomesNoRepMasker/1_4genomesNoRepMasker_config.yaml \
    -t 8 \
    --dry-run
# Expected: Executing: snakemake ... --cores 8  (no --executor, no --jobs)
```

---

### Upload to Anaconda Cloud

```bash
cd /data/toby/EarlGreyParTEA
conda build conda/
anaconda login
anaconda upload /data/toby/miniforge3/conda-bld/noarch/earlgrey-partea-0.1.5-py_0.conda
```

Test the new version in a fresh environment:

```bash
mamba create -n partea_015 -c toby_baril_bio earlgrey-partea=0.1.5
conda activate partea_015
ln -s /data/toby/tools/earlgrey_databases/Libraries/famdb/* /data/toby/miniforge3/envs/partea_015/share/RepeatMasker/Libraries/famdb/
ln -s /data/toby/tools/earlgrey_databases/Libraries/RMRB.embl /data/toby/miniforge3/envs/partea_015/share/RepeatMasker/Libraries/RMRB.embl
ln -s /data/toby/tools/earlgrey_databases/Libraries/RMRBSeqs.embl /data/toby/miniforge3/envs/partea_015/share/RepeatMasker/Libraries/RMRBSeqs.embl
cd /data/toby/miniforge3/envs/partea_015/share/RepeatMasker/
/data/toby/miniforge3/envs/partea_015/bin/perl ./configure
cd /data/toby/testDIR/saturationTests/6_slurmTest
```

## Release v0.1.6 Feature Updates
I would like to add some more features to provide users with useful summaries and outputs. Also, I want to reduce the technical requirements to enable non-specialists to get useful outputs from the pipeline.

Feature list for development:
- Add a `shared_unique_content` rule that plots the number of shared and unique TE families between all species. This can use the cluster file to determine whether a TE family is shared or unique, then plot bars with x-axis as species or assembly, and y-axis as count of shared/unique families, with colour gradients consistent with the plots generated by Earl Grey (i.e. colour patterns used in autoPie.R, with shades for shared vs unique). This rule should also generate the same plots for TE coverage in each genome to show how much of the genome is covered by shared vs unique TEs. A TSV file with the underlying data should also be generated for both family counts and coverage, to allow users to make their own plots if desired.

- Add a `busco_simple_phylo` smk with several rules. First, one that takes as input the genome list, a string argument for an odb lineage set (e.g. `fungi_odb10`), and a string argument for the output prefix. This rule should run BUSCO in genome mode on each input genome using the specified lineage set, then generate a simple summary table with the BUSCO completeness scores for each genome, and a simple bar plot showing the completeness scores across genomes. This will provide users with a quick way to assess the quality of their input genomes and how they compare to each other. In addition, this smk file will also have rules that take BUSCO outputs, finds single copy amino acid sequences with <5% missing data across the whole dataset, makes alignments for each gene, then generates a supermatrix using phykit create_concat, and finally runs fasttree to generate a simple phylogenetic tree based on the BUSCO genes. This will provide users with a quick and easy way to get a phylogenetic tree of their species based on conserved single copy genes, which can be useful for interpreting the TE content in an evolutionary context. If users choose to run this module, the shared_unique_content rule should also generate a plot with assemblies ordered according to the phylogenetic tree, to allow users to see how shared and unique TE content varies across the phylogeny. This should draw the phylogeny as a simple cladogram with the same colour scheme as the shared/unique plots, to visually link the phylogeny with the TE content. Another rule should also show BUSCO scores vs number of TE families discovered in each genome, to allow users to see if there is any relationship between genome quality and TE discovery and to use as a QC diagnostic.

---

## v0.1.6 Implementation Plan

### Overview

Two new optional modules, each activated by a config flag and included in the `Snakefile` only when requested. Both modules ship as dedicated `.smk` files in `rules/` and Python scripts in `scripts/`. Enabling them does **not** change the behaviour of existing rules; all new config keys default to `false`/empty.

| Module | Config flag | Rule file | Required pipeline modes |
|--------|-------------|-----------|------------------------|
| Shared/unique TE content | `run_shared_unique: true` | `rules/shared_unique_content.smk` | `full` only (needs both cluster file and annotation outputs) |
| BUSCO phylogenomics | `run_busco_phylo: true` | `rules/busco_phylo.smk` | Any (uses raw input genomes from config; independent of lib/annotate pipeline) |

---

### New conda dependencies (`conda/meta.yaml`)

The following packages must be added to the `run:` section. They are only exercised when the respective module is enabled, but they must be present in the environment for the `include:` directive to parse without errors.

| Package | Version | Used by | Source channel |
|---------|---------|---------|----------------|
| `busco` | `>=5.4` | `busco_phylo` | `bioconda` |
| `mafft` | `>=7.450` | `busco_phylo` | `bioconda` |
| `phykit` | `>=1.11` | `busco_phylo` | `bioconda` |
| `clipkit` | `>=2.0` | `busco_phylo` | `bioconda` |
| `fasttree` | `>=2.1.11` | `busco_phylo` | `bioconda` |
| `biopython` | `>=1.81` | `busco_phylo` (tree parsing) | `conda-forge` |
| `pandas` | `>=1.5` | `shared_unique_content`, `busco_phylo` | `conda-forge` |

`pandas` is likely already a transitive dependency, but should be listed explicitly. `numpy` and `matplotlib` are already in the recipe.

---

### New config keys

Add the following block to `config/config.yaml` and all three `--generate-config` template strings in `earlGreyParTEA`, `earlGreyParTEA_LibConstruct`, `earlGreyParTEA_AnnotationOnly`, and `scripts/generate_config.py`:

```yaml
# Optional modules (disabled by default)
run_shared_unique: false    # Generate shared/unique TE content plots and tables
                            # Requires pipeline_mode: "full" (needs cluster file + annotations)

run_busco_phylo: false      # Run BUSCO-based phylogenomics module
busco_lineage: ""           # BUSCO ODB lineage set (e.g. "fungi_odb10", "insecta_odb10")
                            # Required when run_busco_phylo: true
busco_prefix: "busco"       # Output prefix for supermatrix and tree files
busco_min_occupancy: 0.95   # Minimum fraction of genomes a gene must be present in
                            # to be included in the supermatrix (default: 0.95 = <5% missing)
```

**Config validation additions in `scripts/on_start_functions.py`:**
- Add `run_shared_unique`, `run_busco_phylo`, `busco_lineage`, `busco_prefix`, `busco_min_occupancy` to the `defaults` dict.
- Emit `[ERROR]` if `run_shared_unique: true` and `pipeline_mode` is `"libconstruct"` (no annotation BED files will exist — cannot determine which genomes contain each TE).
- Emit `[WARNING]` if `run_shared_unique: true` and `pipeline_mode` is `"annotate"` (presence/absence mode will be used; sequence-level divergence not accounted for; cluster-based analysis is more accurate — consider running in `full` mode).
- Emit `[INFO]` if `run_shared_unique: true` and `pipeline_mode` is `"full"` (cluster-based — most accurate).
- Emit `[ERROR]` if `run_busco_phylo: true` and `busco_lineage` is empty.
- Emit `[INFO]` messages describing what each enabled module will do.

---

### Snakefile integration

In `Snakefile`, after the existing `include:` blocks:

```python
RUN_SHARED_UNIQUE = config.get("run_shared_unique", False)
RUN_BUSCO_PHYLO   = config.get("run_busco_phylo",   False)

if RUN_SHARED_UNIQUE:
    include: "rules/shared_unique_content.smk"

if RUN_BUSCO_PHYLO:
    include: "rules/busco_phylo.smk"
```

In the `rule all` for `full` mode, add the new outputs conditionally:

```python
# Inside rule all for "full" mode:
f"{OUTDIR}/combinedLibraries/shared_unique_families.pdf"      if RUN_SHARED_UNIQUE else [],
f"{OUTDIR}/combinedLibraries/shared_unique_families.tsv"      if RUN_SHARED_UNIQUE else [],
f"{OUTDIR}/combinedLibraries/shared_unique_coverage.pdf"      if RUN_SHARED_UNIQUE else [],
f"{OUTDIR}/combinedLibraries/shared_unique_coverage.tsv"      if RUN_SHARED_UNIQUE else [],
f"{OUTDIR}/busco_phylo/{BUSCO_PREFIX}.tree"                   if RUN_BUSCO_PHYLO else [],
f"{OUTDIR}/busco_phylo/{BUSCO_PREFIX}_busco_summary.pdf"      if RUN_BUSCO_PHYLO else [],
f"{OUTDIR}/busco_phylo/{BUSCO_PREFIX}_busco_summary.tsv"      if RUN_BUSCO_PHYLO else [],
# phylo-ordered shared/unique plot only produced when both modules run:
f"{OUTDIR}/combinedLibraries/shared_unique_families_phylo.pdf" if (RUN_SHARED_UNIQUE and RUN_BUSCO_PHYLO) else [],
f"{OUTDIR}/combinedLibraries/shared_unique_coverage_phylo.pdf" if (RUN_SHARED_UNIQUE and RUN_BUSCO_PHYLO) else [],
```

Two global variables read from config and used by the `.smk` files:
```python
BUSCO_LINEAGE = config.get("busco_lineage", "")
BUSCO_PREFIX  = config.get("busco_prefix", "busco")
BUSCO_MIN_OCC = config.get("busco_min_occupancy", 0.95)
```

---

### Module 1: `rules/shared_unique_content.smk`

#### Conceptual approach

Shared/unique classification depends on what information is available, which differs between pipeline modes. The module supports **two detection strategies** selected automatically from `pipeline_mode`:

---

**Strategy A — cluster-based** (`full` mode only)

The cd-hit-est `.clstr` file encodes cross-genome cluster membership. Each sequence in the combined library has a `{species}_` prefix (or `REPMASKER_`/`CUSTOM_` for third-party libraries), following the convention established in `clustering.smk`. Parsing this file gives, for each cluster:

- the set of species that contributed at least one sequence
- the cluster representative (longest sequence, marked with `*` in the file)

*Family classification:*
- A cluster is **unique** to species X if every sequence in it carries the `{species_X}_` prefix.
- A cluster is **shared** if sequences from ≥2 species are present.

*Coverage classification:*
The per-species `.filteredRepeats.gff` (in `{species}_summaryFiles/`) is parsed for coverage. Each GFF row carries a `NAME=` attribute in column 9 (e.g., `NAME=GENOME1_RND-1_FAMILY-34`) — this is the canonical TE family name, consistent across all annotation instances of that family. **The BED file column 4 contains only the TE class/subclass (e.g., `LINE/R2-Hero`) not the family name, so the GFF must be used.** The family name is looked up in `rep_to_cluster` to determine shared vs unique status. Coverage bp = `end − start + 1` (GFF coordinates are 1-based, end-inclusive).

---

**Strategy B — presence/absence-based** (`annotate` mode)

In `annotate` mode the user supplies their own library; no `.clstr` file is produced. Instead, shared/unique status is determined from the annotation outputs alone: for a given named TE family, which genomes have at least one annotated hit in their `.filteredRepeats.gff`?

The `NAME=` attribute in GFF column 9 provides the canonical TE family name. **The BED file must not be used for this — its column 4 is the TE class, not the family name.**

*Family classification:*
- A TE family `NAME=` value is **unique** if it appears in only one species' GFF.
- It is **shared** if it appears in ≥2 species' GFFs.

*Coverage classification:*
Same GFF iteration as Strategy A, but using the presence/absence family→species map instead of a cluster membership dict. Coverage bp = `end − start + 1`.

Important caveat to add as a warning in the TSV and log: presence/absence classification does **not** account for sequence divergence between species. A TE that is genuinely shared but under a different name in each species' annotation will appear as two unique families. Strategy A (cluster-based) is therefore more accurate when available — users should prefer `full` mode when they want shared/unique analysis. Emit `[WARNING]` via `on_start_functions.py` when `run_shared_unique: true` and `pipeline_mode: annotate`.

---

Genome sizes (for computing percentage coverage) can be computed from the `.prep` FASTA files that are always produced during the pipeline in all modes.

#### Implementation

All logic lives in `scripts/shared_unique_plot.py` (analogous to `saturation_plot.py`). It is called from a Snakemake `script:` directive, receiving params from the rule. The rule definition differs between `full` and `annotate` modes; the Snakefile uses `PIPELINE_MODE` to select which rule variant to include.

**Rule: `shared_unique_content`** (`full` mode — cluster-based)

```
Input:
  clstr  — {OUTDIR}/combinedLibraries/combined_all_species.clstrd.fa.clstr
  gffs   — expand("{outdir}/{species}_EarlGrey/{species}_summaryFiles/{species}.filteredRepeats.gff",
                  outdir=OUTDIR, species=SPECIES_LIST)
  preps  — expand("{outdir}/{species}_EarlGrey/{species}.prep",
                  outdir=OUTDIR, species=SPECIES_LIST)

Output:
  fam_pdf   — {OUTDIR}/combinedLibraries/shared_unique_families.pdf
  fam_tsv   — {OUTDIR}/combinedLibraries/shared_unique_families.tsv
  cov_pdf   — {OUTDIR}/combinedLibraries/shared_unique_coverage.pdf
  cov_tsv   — {OUTDIR}/combinedLibraries/shared_unique_coverage.tsv

threads: 1
resources: mem_mb=4000, runtime=30
params:
  species        = SPECIES_LIST
  outdir         = OUTDIR
  detection_mode = "cluster"  # tells the script which code path to use
  has_phylo_tree = False       # always False in this rule; phylo variant is separate
script: "../scripts/shared_unique_plot.py"
```

**Rule: `shared_unique_content_pa`** (`annotate` mode — presence/absence)

Same outputs, different inputs — no `.clstr` file. The per-species GFF files are the sole source of both family names and coverage coordinates.

```
Input:
  gffs   — expand("{outdir}/{species}_EarlGrey/{species}_summaryFiles/{species}.filteredRepeats.gff",
                  outdir=OUTDIR, species=SPECIES_LIST)
  preps  — expand("{outdir}/{species}_EarlGrey/{species}.prep",
                  outdir=OUTDIR, species=SPECIES_LIST)

Output:  (identical to cluster-based rule)

threads: 1
resources: mem_mb=4000, runtime=30
params:
  species        = SPECIES_LIST
  outdir         = OUTDIR
  detection_mode = "presence_absence"  # triggers alternate script code path
  has_phylo_tree = False
script: "../scripts/shared_unique_plot.py"
```

The Snakefile selects the correct rule variant:

```python
if RUN_SHARED_UNIQUE:
    if PIPELINE_MODE == "full":
        include: "rules/shared_unique_content.smk"    # cluster-based variant
    elif PIPELINE_MODE == "annotate":
        include: "rules/shared_unique_content_pa.smk"  # presence/absence variant
    # libconstruct → error already raised in validate_parameters (no annotation BEDs exist)
```

Alternatively, both variants can live in a single `rules/shared_unique_content.smk` with the rule name chosen at parse time; either approach works. The separate-file approach is cleaner for DAG readability.

**Rule: `shared_unique_content_phylo`** (only included when both modules enabled)

Identical to the mode-appropriate non-phylo rule but with the tree as an additional input and `has_phylo_tree=True`. Works in both `full` and `annotate` modes (whichever base rule is active).

```
Input:  (mode-appropriate inputs above) + tree={OUTDIR}/busco_phylo/{BUSCO_PREFIX}.tree
Output:
  fam_phylo_pdf — {OUTDIR}/combinedLibraries/shared_unique_families_phylo.pdf
  cov_phylo_pdf — {OUTDIR}/combinedLibraries/shared_unique_coverage_phylo.pdf
script: "../scripts/shared_unique_plot.py"  # same script, tree_path param set
```

#### `scripts/shared_unique_plot.py` — functional description

The script branches early on `snakemake.params.detection_mode` (`"cluster"` or `"presence_absence"`).

---

**Mode A: `"cluster"` (full mode)**

*Step 1 — parse cluster file*
Re-use `parse_clstr()` from `scripts/cluster_utils.py`. Build:
- `cluster_species: dict[cluster_id → frozenset[str]]` — species set per cluster
- `rep_to_cluster: dict[rep_name → cluster_id]` — map representative name → cluster id

*Step 2 — family counts*
For each species, count:
- `n_unique` = clusters where `cluster_species[cid] == {species}`
- `n_shared` = clusters where `species in cluster_species[cid]` and `len(cluster_species[cid]) > 1`

*Step 4 — coverage*
For each species, iterate through its `.filteredRepeats.gff`. Parse the `NAME=` value from column 9 attributes to get the canonical family name. Look it up in `rep_to_cluster` (the `NAME=` value corresponds directly to the sequence name used when the combined library was built, e.g., `GENOME1_RND-1_FAMILY-34`). Route `end − start + 1` bp (GFF is 1-based, end-inclusive) into `shared_bp` or `unique_bp`.

---

**Mode B: `"presence_absence"` (annotate mode)**

*Step 1 — build family→species map from GFF files*
```python
family_species: dict[str, set[str]]  # family_name → {species, ...}
```
For each species' `.filteredRepeats.gff`, parse the `NAME=` attribute from column 9 to extract the canonical TE family name. Record that species as having observed this family. **Do not use column 4 — that is the TE class/subclass, not the family name.**

*Step 2 — family counts*
For each species:
- `n_unique` = families in `family_species` where `family_species[fam] == {species}`
- `n_shared` = families where `species in family_species[fam]` and `len(family_species[fam]) > 1`

*Step 4 — coverage*
For each GFF row, parse `NAME=` from column 9 to get the family name. Look it up in `family_species`. If it was seen in only this species → `unique_bp += end − start + 1`; if seen in ≥2 species → `shared_bp += end − start + 1`.

A `method` column is added to the TSV output to indicate which detection strategy was used. For Mode B, also write a header comment `# NOTE: presence/absence classification — sequence-level divergence not accounted for` at the top of the TSV.

---

**Shared steps (both modes)**

*Step 3 — genome sizes*
For each `.prep` FASTA, count total bp (sum of sequence lengths). Store as `dict[species → int]`.

*Write outputs*

`shared_unique_families.tsv`:
```
species  shared_families  unique_families  total_families  method
```

`shared_unique_coverage.tsv`:
```
species  shared_bp  unique_bp  total_bp  genome_size_bp  shared_pct  unique_pct  method
```

**Step 5 — plots**

Colour scheme: Two-tone approach matching Earl Grey's aesthetic. Earl Grey's main summary pie chart uses a named colour palette for TE classes via `autoPie.R`. For the shared/unique summary (which collapses across classes), use:
- Shared TEs: `#4477AA` (blue, Earl Grey's DNA element colour, serves as the "shared" indicator)
- Unique TEs: `#BBCCEE` (light blue, same hue but desaturated = unique)
  
Rationale: We cannot know the per-class breakdown at the cross-species level without significantly more parsing complexity. Colour-by-class would require splitting coverage by class for each shared/unique category — this can be added as a v2 enhancement but is out of scope for v0.1.6. Note a comment in the code that class-level breakdown is a future extension.

Both bar plots: horizontal stacked bars, one bar per species, species on y-axis. Shared (darker) stacked below unique (lighter). Add value labels to bars.

**Step 6 — phylo variant** (when `tree_path` is set)
- Read the newick file using `Bio.Phylo` from `biopython`.
- Extract leaf order (tip-to-tip order in the default ladderized layout). This defines the y-axis species order.
- Use `matplotlib` to draw the cladogram as a set of horizontal and vertical line segments in a left-side panel (approx 25% of figure width), with the stacked bars occupying the remaining 75%.
- Scale branch lengths to unit (cladogram style: all tips at same x position). Draw branches in the same dark colour used for the shared bars.
- Species labels appear only on the bar panel (y-tick labels), not on the cladogram, to avoid duplication.

---

### Module 2: `rules/busco_phylo.smk`

All output files go under `{OUTDIR}/busco_phylo/`. Variables `BUSCO_LINEAGE`, `BUSCO_PREFIX`, `BUSCO_MIN_OCC`, `SPECIES_LIST`, `GENOME`, and `OUTDIR` are set at the top of the file (read from `config`), consistent with the pattern in `annotate.smk`.

#### Rule 1: `run_busco`

Runs BUSCO in genome mode per species on the raw input genome (from `config["genome"]`, not the `.prep` file — BUSCO is a genome QC tool and should see the unmodified assembly).

```
Input:  genome = lambda wildcards: GENOME[wildcards.species]
Output: busco_dir = directory("{OUTDIR}/busco_phylo/busco_runs/{species}_{BUSCO_LINEAGE}_busco")
        busco_done = "{OUTDIR}/busco_phylo/busco_runs/{species}_{BUSCO_LINEAGE}_busco/short_summary.specific.{BUSCO_LINEAGE}.{species}_{BUSCO_LINEAGE}_busco.txt"
        # The short_summary file is a reliable sentinel for job completion.
threads: slurm-aware lambda (same pattern as other rules)
resources: mem_mb = lambda wildcards, attempt: 16000 * attempt, runtime = 10080  # 7 days
params:
  lineage = BUSCO_LINEAGE
  outdir  = "{OUTDIR}/busco_phylo/busco_runs"
  prefix  = "{species}_{BUSCO_LINEAGE}_busco"
shell:
  busco -i {input.genome} -m genome -l {params.lineage} \
        -c {threads} -o {params.prefix} --out_path {params.outdir} \
        --offline  # if offline databases are available; otherwise omit
```

Design note: BUSCO downloads lineage data on first run unless `--offline` is specified with pre-downloaded databases. Users with internet access can omit `--offline`. Consider adding a `busco_offline: false` config key (default false = allow download) to give users control. Add this to the config block.

#### Rule 2: `busco_summary_table`

Aggregates all per-species `short_summary.specific.*.txt` files into a summary table and bar chart.

```
Input:  summaries = expand("{OUTDIR}/busco_phylo/busco_runs/{species}_{BUSCO_LINEAGE}_busco/short_summary.specific.{BUSCO_LINEAGE}.{species}_{BUSCO_LINEAGE}_busco.txt", ...)
Output:
  tsv = "{OUTDIR}/busco_phylo/{BUSCO_PREFIX}_busco_summary.tsv"
  pdf = "{OUTDIR}/busco_phylo/{BUSCO_PREFIX}_busco_summary.pdf"
threads: 1
resources: mem_mb=4000, runtime=30
script: "../scripts/busco_summary_plot.py"
```

**`scripts/busco_summary_plot.py`:**
- Parse each `short_summary.specific.*.txt` to extract: species name, Complete%, Complete Single-Copy%, Complete Duplicated%, Fragmented%, Missing%.
- Write TSV with one row per species.
- Generate a horizontally stacked bar chart (standard BUSCO plot style): for each species one bar split into Complete(S), Complete(D), Fragmented, Missing. Use BUSCO's conventional green/yellow/red colour scheme.

#### Checkpoint: `extract_busco_aa`

This is the dynamic step — the number of qualifying BUSCO gene IDs is not known until BUSCO runs complete. A Snakemake `checkpoint` is used here so downstream alignment rules can be dynamically expanded.

```python
checkpoint extract_busco_aa:
    input:
        full_tables = expand(
            "{OUTDIR}/busco_phylo/busco_runs/{species}_{BUSCO_LINEAGE}_busco/run_{BUSCO_LINEAGE}/full_table.tsv",
            OUTDIR=OUTDIR, species=SPECIES_LIST, BUSCO_LINEAGE=BUSCO_LINEAGE
        )
    output:
        aa_dir = directory(f"{OUTDIR}/busco_phylo/aa_sequences"),
        gene_list = f"{OUTDIR}/busco_phylo/filtered_busco_ids.txt"
    threads: 1
    resources: mem_mb=4000, runtime=30
    params:
        species = SPECIES_LIST,
        lineage = BUSCO_LINEAGE,
        min_occupancy = BUSCO_MIN_OCC,
        busco_run_dir = f"{OUTDIR}/busco_phylo/busco_runs"
    script: "../scripts/extract_busco_aa.py"
```

**`scripts/extract_busco_aa.py`:**

Step 1 — find complete single-copy genes per species:
For each species, parse `full_table.tsv`. Keep only rows where column 2 == `"Complete"`. Collect the BUSCO gene IDs.

Step 2 — filter by occupancy:
Count across all species how many have a given gene as "Complete". Keep genes present in `>= ceil(n_species * min_occupancy)` species.
Write `filtered_busco_ids.txt` (one gene ID per line).

Step 3 — extract AA sequences:
For each species and each qualifying gene, locate the `.faa` file at:
`{busco_run_dir}/{species}_{lineage}_busco/run_{lineage}/busco_sequences/single_copy_busco_sequences/{gene_id}.faa`

Rename the FASTA header to just `>{species}` (strip everything else — this is critical for the supermatrix to have consistent taxon labels). Collect all species' sequences for a given gene into a single per-gene file:
`{aa_dir}/{gene_id}_aa.fasta`

This creates one multi-species FASTA per gene ID.

#### Rule 3: `align_busco_gene`

One rule instance per gene (dynamic expansion via checkpoint resolver function).

```
Input:  fasta = f"{OUTDIR}/busco_phylo/aa_sequences/{{gene_id}}_aa.fasta"
Output: aln   = f"{OUTDIR}/busco_phylo/alignments/{{gene_id}}.aln"
threads: 2  # mafft --auto scales modestly; 2 threads is a good default for per-gene alignment
resources: mem_mb=4000, runtime=120
shell:
  mafft --auto --thread {threads} {input.fasta} > {output.aln}
```

#### Checkpoint resolver function

```python
def get_aligned_genes(wildcards):
    # Force Snakemake to wait for the checkpoint output
    checkpoints.extract_busco_aa.get()
    gene_list_path = f"{OUTDIR}/busco_phylo/filtered_busco_ids.txt"
    with open(gene_list_path) as fh:
        genes = [line.strip() for line in fh if line.strip()]
    return expand(f"{OUTDIR}/busco_phylo/alignments/{{gene_id}}.aln", gene_id=genes)
```

#### Rule 4: `fix_aln_headers`

After `mafft`, headers in each `.aln` file still contain the full original BUSCO format (`>speciesName|BUSCO_id:...`). We stripped the BUSCO part during extraction (Step 3 above), so headers should already be clean (`>species`). However, mafft can sometimes modify headers. A quick `sed` pass to keep only the first word after `>` is applied as a safety measure. This can be folded into the `align_busco_gene` shell command rather than a separate rule:

```bash
mafft --auto --thread {threads} {input.fasta} | \
  sed '/^>/s/ .*//g' > {output.aln}
```

Removing the need for a separate `fix_aln_headers` rule.

#### Rule 5: `create_supermatrix`

```
Input:  alignments = get_aligned_genes  (checkpoint resolver)
Output:
  concat_fa    = f"{OUTDIR}/busco_phylo/{BUSCO_PREFIX}_supermatrix.fa"
  concat_parts = f"{OUTDIR}/busco_phylo/{BUSCO_PREFIX}_supermatrix.partition"
  trimmed_fa   = f"{OUTDIR}/busco_phylo/{BUSCO_PREFIX}_supermatrix.fa.clipkit"
threads: 4
resources: mem_mb=16000, runtime=240
params:
  outdir = f"{OUTDIR}/busco_phylo",
  prefix = BUSCO_PREFIX
run:
  import os
  aln_list = f"{params.outdir}/{params.prefix}_alignment_list.txt"
  with open(aln_list, 'w') as fh:
      for aln in sorted(input.alignments):
          fh.write(os.path.abspath(aln) + "\n")
  shell(f"phykit create_concat -a {aln_list} -p {params.outdir}/{params.prefix}_supermatrix")
  shell(f"clipkit {params.outdir}/{params.prefix}_supermatrix.fa -m kpic")
  # clipkit output is automatically named {input}.clipkit
```

Note: `phykit create_concat` requires absolute paths in the alignment list. The `run:` block uses `os.path.abspath()` on each alignment path.

#### Rule 6: `run_fasttree`

```
Input:  trimmed = f"{OUTDIR}/busco_phylo/{BUSCO_PREFIX}_supermatrix.fa.clipkit"
Output: tree    = f"{OUTDIR}/busco_phylo/{BUSCO_PREFIX}.tree"
threads: lambda wildcards: max(1, min(workflow.cores, 8)) if config.get("slurm_mode", False) else max(1, min(workflow.cores // len(SPECIES_LIST), 8))
resources: mem_mb=lambda wildcards, attempt: 16000 * attempt, runtime=1440  # 24 h
shell:
  FastTree -log {output.tree}.log {input.trimmed} > {output.tree}
  # FastTree outputs newick to stdout; log file for diagnostics
```

Note: FastTree uses OpenMP threading automatically. Ensure the conda package is the OpenMP-enabled `fasttree` (the bioconda package is). The number of threads is controlled by the `OMP_NUM_THREADS` environment variable, which Snakemake sets via `{threads}`. Add `export OMP_NUM_THREADS={threads}` before the FastTree call to be explicit.

#### Rule 7: `busco_vs_te_qc`

Scatter plot of BUSCO completeness vs TE families discovered, used as a QC diagnostic.

```
Input:
  busco_tsv = f"{OUTDIR}/busco_phylo/{BUSCO_PREFIX}_busco_summary.tsv"
  sat_tsv   = f"{OUTDIR}/combinedLibraries/saturation_data.tsv"  # only if libconstruct/full mode
  clstr     = f"{OUTDIR}/combinedLibraries/combined_all_species.clstrd.fa.clstr"
Output:
  qc_pdf    = f"{OUTDIR}/busco_phylo/{BUSCO_PREFIX}_busco_vs_te_families.pdf"
  qc_tsv    = f"{OUTDIR}/busco_phylo/{BUSCO_PREFIX}_busco_vs_te_families.tsv"
threads: 1
resources: mem_mb=4000, runtime=30
script: "../scripts/busco_te_qc_plot.py"
```

Design note: this rule is only included when `run_busco_phylo: true` AND `pipeline_mode in ("full", "libconstruct")`. For `annotate` mode (where the cluster file and saturation data don't exist), the rule is skipped. Add a conditional `include:` in the Snakefile or an `if` guard on the `rule all` target list.

**`scripts/busco_te_qc_plot.py`:**
- Read BUSCO summary TSV: map species → complete single-copy %.
- Read saturation data TSV OR parse cluster file to count per-species TE families (using `parse_clstr` from `cluster_utils.py`).
- Scatter plot: x-axis = BUSCO complete single-copy %, y-axis = TE families discovered. Each point labelled with species name. Add a horizontal reference line at the mean.
- Write QC TSV:
  `species  busco_complete_pct  te_family_count`

---

### Script utilities / shared code

To avoid duplicating `parse_clstr()` (which already exists in `saturation_plot.py`), refactor it into `scripts/cluster_utils.py`:

```python
# scripts/cluster_utils.py
def parse_clstr(clstr_file, species_list):
    """Shared cluster parsing utility. Used by saturation_plot.py, shared_unique_plot.py, busco_te_qc_plot.py."""
    ...
```

Update `saturation_plot.py` to import from `cluster_utils` rather than defining the function inline. This is a small refactor with low risk.

---

### File summary

**New files to create:**
- `rules/shared_unique_content.smk`
- `rules/busco_phylo.smk`
- `scripts/shared_unique_plot.py`
- `scripts/busco_summary_plot.py`
- `scripts/extract_busco_aa.py`
- `scripts/busco_te_qc_plot.py`
- `scripts/cluster_utils.py`

**Files to modify:**
- `Snakefile` — module include blocks, rule all targets, global variables
- `scripts/on_start_functions.py` — new config keys in `defaults`, new validation messages and errors, update `show_pipeline_mode_visualization` table
- `scripts/saturation_plot.py` — import `parse_clstr` from `cluster_utils` instead of defining locally
- `scripts/generate_config.py` — add new config keys to all three template strings
- `earlGreyParTEA`, `earlGreyParTEA_LibConstruct`, `earlGreyParTEA_AnnotationOnly` — new config block in the inline `generate_example_config()` function
- `config/config.yaml` — add new config keys with comments
- `conda/meta.yaml` — add new run dependencies

---

### Pipeline mode compatibility matrix

| Module | `full` | `libconstruct` | `annotate` |
|--------|--------|----------------|------------|
| `run_shared_unique` | ✓ Cluster-based (most accurate) | ✗ Error (no annotation BED files) | ✓ Presence/absence-based (with warning) |
| `run_busco_phylo` (BUSCO + tree) | ✓ | ✓ | ✓ |
| `busco_vs_te_qc` (sub-rule) | ✓ | ✓ | ✗ Skipped (omitted from rule all) |
| `shared_unique_content_phylo` (combined plot) | ✓ (if both enabled) | ✗ | ✓ (if both enabled, presence/absence method) |

---

### Resource requirements summary (for SLURM mode)

| Rule | threads | mem_mb | runtime (min) |
|------|---------|--------|---------------|
| `run_busco` | `--threads` value | 16 000 × attempt | 10 080 (7 days) |
| `busco_summary_table` | 1 | 4 000 | 30 |
| `extract_busco_aa` (checkpoint) | 1 | 4 000 | 60 |
| `align_busco_gene` | 2 | 4 000 | 120 per gene |
| `create_supermatrix` | 4 | 16 000 | 240 |
| `run_fasttree` | up to 8 | 16 000 × attempt | 1 440 (24 h) |
| `busco_vs_te_qc` | 1 | 4 000 | 30 |
| `shared_unique_content` | 1 | 4 000 | 30 |
| `shared_unique_content_phylo` | 1 | 4 000 | 30 |

---

### Testing checklist

#### Unit-level tests (dry run / small dataset)

- [ ] `run_shared_unique: false` (default) — confirm no new rules appear in DAG
- [ ] `run_busco_phylo: false` (default) — confirm no new rules appear in DAG
- [ ] `run_shared_unique: true` in `libconstruct` mode — confirm error exit with clear message
- [ ] `run_shared_unique: true` in `annotate` mode — confirm `[WARNING]` printed and `shared_unique_content_pa` rule appears in DAG (not the cluster-based rule)
- [ ] `run_busco_phylo: true` with empty `busco_lineage` — confirm error exit with clear message
- [ ] `run_shared_unique: true` in `full` mode — dry run shows `shared_unique_content` in DAG, downstream targets correct
- [ ] `run_busco_phylo: true` in `full` mode — dry run shows full BUSCO rule chain in DAG including checkpoint

#### Integration tests (real data)

- [ ] Run `busco_summary_table` on an existing BUSCO output; confirm TSV has correct species × BUSCO score columns
- [ ] Run `extract_busco_aa` checkpoint with the Z. tritici 4-genome test set; confirm `filtered_busco_ids.txt` contains only genes present in ≥95% of genomes; confirm per-gene FASTA files have `>species` headers
- [ ] Run full `busco_phylo` chain on 4-genome test set; confirm `{prefix}.tree` is a valid newick; confirm tree has 4 leaves matching the species names exactly
- [ ] Run `shared_unique_content` (cluster mode) on existing 4-genome annotation outputs; confirm TSV `method` column shows `cluster`; confirm total families ≈ count from saturation data; shared+unique = total; PDFs open without error
- [ ] Run `shared_unique_content_pa` (presence/absence mode) on the same 4-genome outputs by switching to `annotate` mode with the clustered library as `annotation_library`; confirm TSV `method` shows `presence_absence`; confirm `[WARNING]` was printed during config validation; confirm family counts are similar to (but may differ from) cluster mode
- [ ] Enable both modules in `full` mode; confirm `shared_unique_content_phylo` PDF shows a cladogram panel with correct species order derived from the tree
- [ ] Enable both modules in `annotate` mode; confirm `shared_unique_content_phylo` uses presence/absence data with correct phylo ordering
- [ ] `busco_vs_te_qc` scatter plot: confirm one point per species, x-axis within 0–100, y-axis ≥ 0

#### Regression tests

- [ ] Run standard `full` mode without either new module enabled — confirm no change in behaviour, outputs, or runtime
- [ ] Run `libconstruct` and `annotate` modes without new modules — confirm no regression

---

## Test Release v0.1.6

### Version bump checklist

Before building, update the version string in all relevant files:

- [x] `earlGreyParTEA` — `VERSION="0.1.6"`
- [x] `earlGreyParTEA_LibConstruct` — `VERSION="0.1.6"`
- [x] `earlGreyParTEA_AnnotationOnly` — `VERSION="0.1.6"`
- [x] `conda/meta.yaml` — `version: "0.1.6"` and update `sha256` after tagging

### Build and install

```bash
cd /data/toby/EarlGreyParTEA
conda-build purge-all
conda build conda/
conda create -n test_016 --use-local earlgrey-partea
conda activate test_016
```

Configure RepeatMasker in the new environment:

```bash
ln -s /data/toby/tools/earlgrey_databases/Libraries/famdb/* /data/toby/miniforge3/envs/test_016/share/RepeatMasker/Libraries/famdb/
ln -s /data/toby/tools/earlgrey_databases/Libraries/RMRB.embl /data/toby/miniforge3/envs/test_016/share/RepeatMasker/Libraries/RMRB.embl
ln -s /data/toby/tools/earlgrey_databases/Libraries/RMRBSeqs.embl /data/toby/miniforge3/envs/test_016/share/RepeatMasker/Libraries/RMRBSeqs.embl
cd /data/toby/miniforge3/envs/test_016/share/RepeatMasker/
/data/toby/miniforge3/envs/test_016/bin/perl ./configure
cd /data/toby/EarlGreyParTEA
```

---

### Test 0: Unit tests and import smoke test

Run the pytest suite to verify all helper functions behave correctly before touching any
real data:

```bash
cd /data/toby/EarlGreyParTEA
python -m pytest tests/ -v
```

Confirm:
- All tests in `tests/test_cluster_utils.py`, `tests/test_shared_unique_plot.py`, and
  `tests/test_busco_modules.py` pass
- No `ImportError` or `NameError` at collection time

Verify that importing any new script as a plain Python module does **not** execute
`main()` (i.e. does not crash with `NameError: name 'snakemake' is not defined`):

```bash
python -c "
import sys; sys.path.insert(0, 'scripts')
from shared_unique_plot import _classify_te, CLASS_PALETTE, LIGHT_PALETTE, CLASS_ORDER
from busco_summary_plot import _parse_busco_summary
from busco_te_qc_plot import _load_coverage_tsv
from cluster_utils import parse_clstr
import extract_busco_aa
print('ALL IMPORTS CLEAN')
"
```

Confirm:
- Prints `ALL IMPORTS CLEAN` with no errors

Spot-check a few TE classification mappings:

```bash
python -c "
import sys; sys.path.insert(0, 'scripts')
from shared_unique_plot import _classify_te, _lighten, CLASS_PALETTE, LIGHT_PALETTE
cases = [
    ('LTR/Gypsy','LTR'), ('LINE/R2-Hero','LINE'), ('LINE/Penelope','Penelope'),
    ('SINE/tRNA','SINE'), ('DNA/TIR','DNA'), ('RC/Helitron','Rolling Circle'),
    ('Simple_repeat','Other'), ('Unknown','Unclassified'),
]
for te, expected in cases:
    got = _classify_te(te)
    print(f'  {te:25s} -> {got:15s} [{'OK' if got == expected else 'FAIL (expected '+expected+')'}]')
# Confirm lightening produces a different (lighter) colour
for cls in CLASS_PALETTE:
    assert CLASS_PALETTE[cls] != LIGHT_PALETTE[cls], f'Lightening did nothing for {cls}'
print('Lightening: OK')
"
```

Confirm:
- Every row ends with `[OK]`
- `LINE/Penelope` maps to `Penelope`, not `LINE` (Penelope-before-LINE priority)
- `Lightening: OK`

---

### Test 1: Config validation error cases

These tests check `on_start_functions.py` validation without running any real tools.
Use `--dry-run` so Snakemake exits after parsing.

**Error: `run_shared_unique: true` in `libconstruct` mode**

```bash
mkdir -p /data/toby/testDIR/0.1.6tests/err_tests

cat > /data/toby/testDIR/0.1.6tests/err_tests/bad_shared_libconstruct.yaml << 'EOF'
pipeline_mode: libconstruct
output_dir: /data/toby/testDIR/0.1.6tests/err_tests/out
genome:
  IPO323: /data/toby/testDIR/saturationTests/1_4genomesNoRepMasker/IPO323.fa
species:
  - IPO323
run_shared_unique: true
run_busco_phylo: false
busco_lineage: ""
EOF

earlGreyParTEA_LibConstruct \
    -c /data/toby/testDIR/0.1.6tests/err_tests/bad_shared_libconstruct.yaml \
    -t 4 \
    --dry-run
# Expected: [ERROR] run_shared_unique=true requires annotation outputs ... not produced in 'libconstruct' mode
```

**Error: `run_busco_phylo: true` with empty `busco_lineage`**

```bash
cat > /data/toby/testDIR/0.1.6tests/err_tests/bad_busco_no_lineage.yaml << 'EOF'
pipeline_mode: full
output_dir: /data/toby/testDIR/0.1.6tests/err_tests/out
genome:
  IPO323: /data/toby/testDIR/saturationTests/1_4genomesNoRepMasker/IPO323.fa
species:
  - IPO323
run_shared_unique: false
run_busco_phylo: true
busco_lineage: ""
EOF

earlGreyParTEA \
    -c /data/toby/testDIR/0.1.6tests/err_tests/bad_busco_no_lineage.yaml \
    -t 4 \
    --dry-run
# Expected: [ERROR] run_busco_phylo=true requires 'busco_lineage' to be set
```

**Error: `busco_min_occupancy` out of range**

```bash
cat > /data/toby/testDIR/0.1.6tests/err_tests/bad_min_occ.yaml << 'EOF'
pipeline_mode: full
output_dir: /data/toby/testDIR/0.1.6tests/err_tests/out
genome:
  IPO323: /data/toby/testDIR/saturationTests/1_4genomesNoRepMasker/IPO323.fa
species:
  - IPO323
run_shared_unique: false
run_busco_phylo: true
busco_lineage: sordariomycetes_odb10
busco_min_occupancy: 1.5
EOF

earlGreyParTEA \
    -c /data/toby/testDIR/0.1.6tests/err_tests/bad_min_occ.yaml \
    -t 4 \
    --dry-run
# Expected: [ERROR] busco_min_occupancy=1.5 is out of range
```

**Warning: `run_shared_unique: true` in `annotate` mode (should warn, not error)**

```bash
cat > /data/toby/testDIR/0.1.6tests/err_tests/warn_shared_annotate.yaml << 'EOF'
pipeline_mode: annotate
output_dir: /data/toby/testDIR/0.1.6tests/err_tests/out
genome:
  IPO323: /data/toby/testDIR/saturationTests/1_4genomesNoRepMasker/IPO323.fa
  s_1A5: /data/toby/testDIR/saturationTests/1_4genomesNoRepMasker/1A5.fa
species:
  - IPO323
  - s_1A5
annotation_library: /data/toby/EarlGreyParTEA/data/dummy.fasta
run_shared_unique: true
run_busco_phylo: false
busco_lineage: ""
EOF

earlGreyParTEA_AnnotationOnly \
    -c /data/toby/testDIR/0.1.6tests/err_tests/warn_shared_annotate.yaml \
    -t 4 \
    --dry-run
# Expected: [WARNING] presence/absence strategy will be used
# Expected: dry run proceeds (not an error)
# Expected: rule all targets include sharedUniqueContent/ outputs
```

**DAG check: no new rules when both modules disabled**

```bash
earlGreyParTEA \
    -c /data/toby/testDIR/saturationTests/1_4genomesNoRepMasker/1_4genomesNoRepMasker_config.yaml \
    -t 4 \
    --dry-run 2>&1 | grep -E "shared_unique|busco"
# Expected: no output (no shared_unique or busco rules in DAG)
```

---

### Test 2: `shared_unique_content` — cluster-based mode (full pipeline)

This test reuses the completed output from an earlier full pipeline run on the 4
Z. tritici genomes. If a prior full-mode result directory already exists, start
from there to save compute time.

```bash
mkdir -p /data/toby/testDIR/0.1.6tests/2_shared_unique_cluster

# Create config with run_shared_unique enabled, full mode
cat > /data/toby/testDIR/0.1.6tests/2_shared_unique_cluster/config.yaml << 'EOF'
pipeline_mode: full
output_dir: /data/toby/testDIR/saturationTests/6_slurmTest
genome:
  s_1A5: /data/toby/testDIR/saturationTests/1_4genomesNoRepMasker/1A5.fa
  Aus01: /data/toby/testDIR/saturationTests/1_4genomesNoRepMasker/Aus01.fa
  IPO323: /data/toby/testDIR/saturationTests/1_4genomesNoRepMasker/IPO323.fa
  YEQ92: /data/toby/testDIR/saturationTests/1_4genomesNoRepMasker/YEQ92.fa
species:
- s_1A5
- Aus01
- IPO323
- YEQ92
run_shared_unique: true
run_busco_phylo: false
busco_lineage: ""
busco_min_occupancy: 0.95
EOF

# Dry run first to confirm correct rules are in the DAG
earlGreyParTEA \
    -c /data/toby/testDIR/0.1.6tests/2_shared_unique_cluster/config.yaml \
    -t 8 \
    --dry-run 2>&1 | grep -E "shared_unique|rule all"
```

Confirm the dry run shows:
- `shared_unique_plot` in the DAG (cluster-based, not `pa`)
- `shared_unique_plot_phylo` is **not** in the DAG (no tree available)
- Targets include `sharedUniqueContent/shared_unique_families.pdf` and `.tsv`

Run the full pipeline:

```bash
earlGreyParTEA \
    -c /data/toby/testDIR/0.1.6tests/2_shared_unique_cluster/config.yaml \
    -t 8
```

Check outputs:

```bash
OUT=/data/toby/testDIR/saturationTests/6_slurmTest
# Confirm all four expected files exist
ls -lh $OUT/sharedUniqueContent/

# TSV: families
head -1 $OUT/sharedUniqueContent/shared_unique_families.tsv
# Expected header includes: species, shared_families, unique_families, then per-class columns

# TSV: coverage
head -1 $OUT/sharedUniqueContent/shared_unique_coverage.tsv
# Expected header includes: species, genome_size_bp, shared_bp, unique_bp, shared_pct, unique_pct, per-class columns

# Sanity: shared + unique should equal total for each species (within rounding)
python -c "
import csv
with open('$OUT/sharedUniqueContent/shared_unique_families.tsv') as f:
    for row in csv.DictReader(f, delimiter='\t'):
        s = int(row['shared_families']); u = int(row['unique_families'])
        print(f\"{row['species']:10s}  shared={s:4d}  unique={u:4d}  total={s+u:4d}\")
"
```

Confirm:
- All four PDFs and TSVs exist in `sharedUniqueContent/`
- Both PDFs open without error
- Per-species `shared + unique = total` (row-wise consistency)
- Bar colour legend shows TE-class colour swatches (LTR, LINE, SINE, DNA, Rolling Circle, Penelope, Other, Unclassified) and shade indicators for shared vs unique

---

### Test 3: `shared_unique_content` — presence/absence mode (annotate pipeline)

Use the clustered library from Test 2 as the `annotation_library` and switch to
`annotate` mode to test the presence/absence fallback path.

```bash
mkdir -p /data/toby/testDIR/0.1.6tests/3_shared_unique_pa

CLUSTERED_LIB=/data/toby/testDIR/saturationTests/6_slurmTest/combinedLibraries/combined_all_species.clstrd.fa

cat > /data/toby/testDIR/0.1.6tests/3_shared_unique_pa/config.yaml << 'EOF'
pipeline_mode: annotate
output_dir: /data/toby/testDIR/0.1.6tests/3_shared_unique_pa/out
genome:
  s_1A5: /data/toby/testDIR/saturationTests/1_4genomesNoRepMasker/1A5.fa
  Aus01: /data/toby/testDIR/saturationTests/1_4genomesNoRepMasker/Aus01.fa
  IPO323: /data/toby/testDIR/saturationTests/1_4genomesNoRepMasker/IPO323.fa
  YEQ92: /data/toby/testDIR/saturationTests/1_4genomesNoRepMasker/YEQ92.fa
species:
  - IPO323
  - s_1A5
  - YEQ92
  - Aus01
annotation_library: /data/toby/testDIR/saturationTests/6_slurmTest/combinedLibraries/combined_all_species.clstrd.fa
run_shared_unique: true
run_busco_phylo: false
busco_lineage: ""
busco_min_occupancy: 0.95
EOF

earlGreyParTEA_AnnotationOnly \
    -c /data/toby/testDIR/0.1.6tests/3_shared_unique_pa/config.yaml \
    -t 8
```

Confirm:
- `[WARNING]` about presence/absence strategy printed at startup
- Rule `shared_unique_pa_plot` (not `shared_unique_plot`) appears in the log
- `sharedUniqueContent/` contains the same four output files
- TSV `method` column (if present) reads `presence_absence`
- Family counts are broadly similar to Test 2 but may differ (families annotated under
  different names in different species appear as unique)

---

### Test 4: `busco_phylo` — BUSCO completeness only

Test just the BUSCO summary step (stop after `busco_summary_table`) to confirm
BUSCO runs and the completeness plot is produced without running the full
phylogenetics chain yet.

```bash
mkdir -p /data/toby/testDIR/0.1.6tests/4_busco_phylo

cat > /data/toby/testDIR/0.1.6tests/4_busco_phylo/config.yaml << 'EOF'
pipeline_mode: annotate
output_dir: /data/toby/testDIR/0.1.6tests/4_busco_phylo/out
genome:
  IPO323: /data/toby/testDIR/saturationTests/1_4genomesNoRepMasker/IPO323.fa
  s_1A5:  /data/toby/testDIR/saturationTests/1_4genomesNoRepMasker/1A5.fa
  YEQ92:  /data/toby/testDIR/saturationTests/1_4genomesNoRepMasker/YEQ92.fa
  Aus01:  /data/toby/testDIR/saturationTests/1_4genomesNoRepMasker/Aus01.fa
species:
  - IPO323
  - s_1A5
  - YEQ92
  - Aus01
annotation_library: /data/toby/testDIR/0.1.6tests/3_shared_unique_pa/out/combinedLibraries/combined_all_species.clstrd.fa
run_shared_unique: false
run_busco_phylo: true
busco_lineage: fungi_odb10
busco_prefix: busco
busco_min_occupancy: 0.5
EOF

# First, dry run to confirm the full rule chain is in the DAG
earlGreyParTEA_AnnotationOnly \
    -c /data/toby/testDIR/0.1.6tests/4_busco_phylo/config.yaml \
    -t 8 \
    --dry-run 2>&1 | grep -E "run_busco|busco_summary|extract_busco|align_busco|create_super|run_fasttree|busco_te_qc"
```

Confirm the dry run shows all 7 rules:
- `run_busco` (one per species)
- `busco_summary_table`
- `extract_busco_aa` (checkpoint)
- `align_busco_gene` (will be resolved after checkpoint)
- `create_supermatrix`
- `run_fasttree`
- `busco_te_qc` **not** in DAG — this requires `shared_unique_coverage.tsv`, which is
  only produced when `run_shared_unique` is also enabled (see Test 5)

Run only through the summary step first to verify BUSCO works:

```bash
earlGreyParTEA_AnnotationOnly \
    -c /data/toby/testDIR/0.1.6tests/4_busco_phylo/config.yaml \
    -t 8 \
    --snakemake-args "--until busco_summary_table"
```

Check BUSCO summary outputs:

```bash
OUT=/data/toby/testDIR/0.1.6tests/4_busco_phylo/out

# Confirm short_summary files exist for each species
ls $OUT/IPO323_EarlGrey/IPO323_busco/short_summary.specific.fungi_odb10.IPO323_busco.txt
ls $OUT/s_1A5_EarlGrey/s_1A5_busco/short_summary.specific.fungi_odb10.s_1A5_busco.txt

# Confirm completeness TSV and PDF
ls -lh $OUT/buscoPhylo/busco_completeness.pdf
cat $OUT/buscoPhylo/busco_completeness.tsv
```

Confirm:
- All 4 short_summary files exist
- `busco_completeness.tsv` has one row per species with `single`, `duplicated`,
  `fragmented`, `missing`, `total` columns
- `single + duplicated + fragmented + missing == total` for every row (within 1 due
  to rounding)
- `busco_completeness.pdf` opens and shows 4 horizontal stacked bars

---

### Test 5: `busco_phylo` — full phylogenomics chain

Continue the run in Test 4 to completion to produce the supermatrix and tree.

```bash
earlGreyParTEA_AnnotationOnly \
    -c /data/toby/testDIR/0.1.6tests/4_busco_phylo/config.yaml \
    -t 8
```

Check phylogenomics outputs:

```bash
OUT=/data/toby/testDIR/0.1.6tests/4_busco_phylo/out

# Occupancy table: gene × species presence/absence
head $OUT/buscoPhylo/busco_gene_occupancy.tsv

# Per-gene single-copy FASTAs (confirm multiple .faa files)
ls $OUT/buscoPhylo/busco_genes/ | head -5
echo "Gene count: $(ls $OUT/buscoPhylo/busco_genes/*.faa | wc -l)"

# Trimmed alignments
ls $OUT/buscoPhylo/aligned/ | head -5

# Supermatrix
grep "^>" $OUT/buscoPhylo/supermatrix.fa
# Expected: 4 lines (one per species), names matching species list exactly

# Check species_tree.nwk is a non-empty newick
cat $OUT/buscoPhylo/species_tree.nwk
python -c "
from Bio import Phylo; import io
with open('$OUT/buscoPhylo/species_tree.nwk') as f: nwk = f.read()
tree = Phylo.read(io.StringIO(nwk), 'newick')
leaves = [c.name for c in tree.get_terminals()]
print('Leaves:', leaves)
assert len(leaves) == 4, f'Expected 4 leaves, got {len(leaves)}'
print('Tree OK')
"
```

Confirm:
- `busco_gene_occupancy.tsv` lists genes with occupancy ≥ 0.5 (column for each species)
- At least one gene FASTA exists in `busco_genes/` (realistic for sordariomycetes)
- Supermatrix has exactly 4 `>species` sequences
- `species_tree.nwk` parses as a valid newick with 4 leaves matching species names

---

### Test 6: Both modules enabled — phylo-ordered shared/unique plots and BUSCO vs TE QC

This is the combined integration test. Both `run_shared_unique` and `run_busco_phylo`
must be `true` simultaneously. The tree produced by Test 5 is reused via Snakemake's
dependency tracking.

```bash
mkdir -p /data/toby/testDIR/0.1.6tests/6_combined

cat > /data/toby/testDIR/0.1.6tests/6_combined/config.yaml << 'EOF'
pipeline_mode: full
output_dir: /data/toby/testDIR/0.1.6tests/6_combined/out
genome:
  IPO323: /data/toby/testDIR/saturationTests/1_4genomesNoRepMasker/IPO323.fa
  s_1A5:  /data/toby/testDIR/saturationTests/1_4genomesNoRepMasker/1A5.fa
  YEQ92:  /data/toby/testDIR/saturationTests/1_4genomesNoRepMasker/YEQ92.fa
  Aus01:  /data/toby/testDIR/saturationTests/1_4genomesNoRepMasker/Aus01.fa
species:
  - IPO323
  - s_1A5
  - YEQ92
  - Aus01
run_shared_unique: true
run_busco_phylo: true
busco_lineage: sordariomycetes_odb10
busco_prefix: busco
busco_min_occupancy: 0.5
EOF

# Dry run: confirm phylo-ordered rules and busco_te_qc appear
earlGreyParTEA \
    -c /data/toby/testDIR/0.1.6tests/6_combined/config.yaml \
    -t 8 \
    --dry-run 2>&1 | grep -E "shared_unique_plot_phylo|shared_unique_pa_plot_phylo|busco_te_qc"
```

Confirm the dry run shows:
- `shared_unique_plot_phylo` in the DAG (not the non-phylo version)
- `busco_te_qc` in the DAG (requires both `species_tree.nwk` and `shared_unique_coverage.tsv`)

Run the full pipeline:

```bash
earlGreyParTEA \
    -c /data/toby/testDIR/0.1.6tests/6_combined/config.yaml \
    -t 8
```

Check all expected outputs:

```bash
OUT=/data/toby/testDIR/0.1.6tests/6_combined/out

echo "=== sharedUniqueContent ==="
ls -lh $OUT/sharedUniqueContent/
# Expected: shared_unique_families.pdf, .tsv, shared_unique_coverage.pdf, .tsv
#           shared_unique_families_phylo.pdf, shared_unique_coverage_phylo.pdf

echo "=== buscoPhylo ==="
ls -lh $OUT/buscoPhylo/
# Expected: busco_completeness.pdf, .tsv, busco_genes/, aligned/,
#           supermatrix.fa, species_tree.nwk, busco_te_qc.pdf, .tsv

# Verify busco_te_qc TSV has one row per species
cat $OUT/buscoPhylo/busco_te_qc.tsv
# Expected: species, busco_completeness_pct, total_te_pct, dominant_class, genome_size_bp columns
```

Confirm:
- `sharedUniqueContent/` contains 6 files (4 plain + 2 phylo PDFs)
- Phylo PDFs show a cladogram panel to the left of the stacked bars, species ordered
  top-to-bottom matching tree leaf order
- `busco_te_qc.pdf` scatter has one point per species; x-axis 0–100 (BUSCO %),
  y-axis ≥ 0 (total TE %); points coloured by dominant TE class using the Earl Grey palette
- `busco_te_qc.tsv` has exactly 4 data rows

---

### Test 7: Regression — no new modules (confirm no behaviour change)

Re-run the standard 4-genome full pipeline without any new modules enabled, using an
output directory separate from the new-module tests, and confirm outputs are identical
to a v0.1.5 run.

```bash
mkdir -p /data/toby/testDIR/0.1.6tests/7_regression

cat > /data/toby/testDIR/0.1.6tests/7_regression/config.yaml << 'EOF'
pipeline_mode: full
output_dir: /data/toby/testDIR/0.1.6tests/7_regression/out
genome:
  IPO323: /data/toby/testDIR/saturationTests/1_4genomesNoRepMasker/IPO323.fa
  s_1A5:  /data/toby/testDIR/saturationTests/1_4genomesNoRepMasker/1A5.fa
  YEQ92:  /data/toby/testDIR/saturationTests/1_4genomesNoRepMasker/YEQ92.fa
  Aus01:  /data/toby/testDIR/saturationTests/1_4genomesNoRepMasker/Aus01.fa
species:
  - IPO323
  - s_1A5
  - YEQ92
  - Aus01
run_shared_unique: false
run_busco_phylo: false
busco_lineage: ""
EOF

earlGreyParTEA \
    -c /data/toby/testDIR/0.1.6tests/7_regression/config.yaml \
    -t 24

# Confirm no buscoPhylo or sharedUniqueContent directories were created
ls /data/toby/testDIR/0.1.6tests/7_regression/out/
# Expected: combinedLibraries/, IPO323_EarlGrey/, s_1A5_EarlGrey/, YEQ92_EarlGrey/,
#           Aus01_EarlGrey/, workflow_visualization/ — NO buscoPhylo/ or sharedUniqueContent/
```

Confirm:
- No `buscoPhylo/` or `sharedUniqueContent/` directories created
- `saturation_plot.pdf` and `saturation_data.tsv` still exist in `combinedLibraries/`
- Per-species `summaryFiles/` outputs all present as expected

### Test 8: Full options with 10 genomes
This is the ultimate test of all new features together on a larger dataset. It will
require a long runtime, so it can be run as a final confirmation before release, but is not essential for the v0.1.6 release itself.

Make conda environment from development branch:
```bash
cd /data/toby/EarlGreyParTEA
conda-build purge-all
mamba env remove -n test_016
conda build conda/
conda create -n test_016 --use-local earlgrey-partea
conda activate test_016

ln -s /data/toby/tools/earlgrey_databases/Libraries/famdb/* /data/toby/miniforge3/envs/test_016/share/RepeatMasker/Libraries/famdb/
ln -s /data/toby/tools/earlgrey_databases/Libraries/RMRB.embl /data/toby/miniforge3/envs/test_016/share/RepeatMasker/Libraries/RMRB.embl
ln -s /data/toby/tools/earlgrey_databases/Libraries/RMRBSeqs.embl /data/toby/miniforge3/envs/test_016/share/RepeatMasker/Libraries/RMRBSeqs.embl
cd /data/toby/miniforge3/envs/test_016/share/RepeatMasker/
/data/toby/miniforge3/envs/test_016/bin/perl ./configure

cd /data/toby/EarlGreyParTEA
```

Make a directory with a few genomes in and make the config:

```bash
mkdir -p /data/toby/testDIR/0.1.6tests/8_full/0_genomes/

cp /legserv/NGS_data/Zymoseptoria/Zt_Reference_genomes/19Pangenome_genomes/Arg00/Arg00.fa /data/toby/testDIR/0.1.6tests/8_full/0_genomes/
cp /legserv/NGS_data/Zymoseptoria/Zt_Reference_genomes/19Pangenome_genomes/Aus01/Aus01.fa /data/toby/testDIR/0.1.6tests/8_full/0_genomes/
cp /legserv/NGS_data/Zymoseptoria/Zt_Reference_genomes/19Pangenome_genomes/CNR93/CNR93.fa /data/toby/testDIR/0.1.6tests/8_full/0_genomes/
cp /legserv/NGS_data/Zymoseptoria/Zt_Reference_genomes/19Pangenome_genomes/I93/I93.fa /data/toby/testDIR/0.1.6tests/8_full/0_genomes/
cp /legserv/NGS_data/Zymoseptoria/Zt_Reference_genomes/19Pangenome_genomes/IPO323/Zymoseptoria_tritici.MG2.dna.toplevel.mt+.fa /data/toby/testDIR/0.1.6tests/8_full/0_genomes/
cp /legserv/NGS_data/Zymoseptoria/Zt_Reference_genomes/19Pangenome_genomes/ISY92/ISY92.fa /data/toby/testDIR/0.1.6tests/8_full/0_genomes/
cp /legserv/NGS_data/Zymoseptoria/Zt_Reference_genomes/19Pangenome_genomes/OregS90/OregS90.fa /data/toby/testDIR/0.1.6tests/8_full/0_genomes/
cp /legserv/NGS_data/Zymoseptoria/Zt_Reference_genomes/19Pangenome_genomes/UR95/UR95.fa /data/toby/testDIR/0.1.6tests/8_full/0_genomes/
cp /legserv/NGS_data/Zymoseptoria/Zt_Reference_genomes/19Pangenome_genomes/YEQ92/YEQ92.fa /data/toby/testDIR/0.1.6tests/8_full/0_genomes/
cp /legserv/NGS_data/Zymoseptoria/Zt_Reference_genomes/19Pangenome_genomes/KE94/KE94.fa /data/toby/testDIR/0.1.6tests/8_full/0_genomes/

earlGreyParTEA \
  --generate-config /data/toby/testDIR/0.1.6tests/8_full/manyGenomesTest.yaml \
  --genome-dir /data/toby/testDIR/0.1.6tests/8_full/0_genomes/ \
  --output-dir /data/toby/testDIR/0.1.6tests/8_full/out/
```

Edit the config to enable lots of things, including slurm!

```yaml
# EarlGrey Pangenome Pipeline Configuration
# Full Mode: Complete library construction and annotation

# Input genomes
genome:
  Arg00: /data/toby/testDIR/0.1.6tests/8_full/0_genomes/Arg00.fa
  Aus01: /data/toby/testDIR/0.1.6tests/8_full/0_genomes/Aus01.fa
  CNR93: /data/toby/testDIR/0.1.6tests/8_full/0_genomes/CNR93.fa
  I93: /data/toby/testDIR/0.1.6tests/8_full/0_genomes/I93.fa
  ISY92: /data/toby/testDIR/0.1.6tests/8_full/0_genomes/ISY92.fa
  KE94: /data/toby/testDIR/0.1.6tests/8_full/0_genomes/KE94.fa
  OregS90: /data/toby/testDIR/0.1.6tests/8_full/0_genomes/OregS90.fa
  UR95: /data/toby/testDIR/0.1.6tests/8_full/0_genomes/UR95.fa
  YEQ92: /data/toby/testDIR/0.1.6tests/8_full/0_genomes/YEQ92.fa
  Zymoseptoria_tritici_MG2_dna_toplevel_mt_: /data/toby/testDIR/0.1.6tests/8_full/0_genomes/Zymoseptoria_tritici.MG2.dna.toplevel.mt+.fa

species:
  - Arg00
  - Aus01
  - CNR93
  - I93
  - ISY92
  - KE94
  - OregS90
  - UR95
  - YEQ92
  - Zymoseptoria_tritici_MG2_dna_toplevel_mt_

# Output directory for all results
output_dir: /data/toby/testDIR/0.1.6tests/8_full/out

# Pipeline mode (do not change for earlGreyParTEA)
pipeline_mode: "full"

# Initial masking with known repeats (choose ONE or leave both empty)
repeatmasker_species: ""  # e.g., "fungi", "arthropoda", "viridiplantae"
custom_library: ""        # path to custom TE library in fasta format

# Library construction parameters
iterations: 10           # Number of BLAST-extend-align cycles
flank: 1000             # Flanking basepairs to extract
max_consensus_seqs: 20  # Max sequences for consensus building
min_consensus_seqs: 3   # Min sequences required for consensus

# Clustering options (for combining TE libraries from multiple genomes)
skip_clustering: false          # Set to true to skip clustering (just concatenate)
clustering_identity: 0.8        # cd-hit sequence identity threshold (0.0-1.0)
clustering_coverage: 0.8        # cd-hit alignment coverage (0.0-1.0)

# Output options
softmask: true         # Generate softmasked genome for each input
margin: false           # Remove short TE sequences (<100bp)
run_heliano: true       # Run HELIANO for Helitron detection

# Workflow visualization (requires graphviz/dot installed)
generate_dag: true      # Generate workflow DAG visualizations
dag_format: "svg"       # Options: "svg", "png", "pdf"

# Saturation analysis options
saturation_permutations: 100  # Number of random genome-addition permutations to average
                              # over for the TE family saturation plot. Higher values
                              # give smoother confidence intervals but increase runtime.

# ---------------------------------------------------------------------------
# Optional analysis modules (v0.1.6+)
# ---------------------------------------------------------------------------

# Shared/unique TE content analysis (full and annotate modes only)
run_shared_unique: true

# BUSCO-based phylogenomics
run_busco_phylo: true
busco_lineage: "fungi_odb10"        # REQUIRED if run_busco_phylo: true  (e.g. "fungi_odb10")
busco_prefix: "busco"   # Prefix for BUSCO run directory names
busco_min_occupancy: 0.1  # Min fraction of species a gene must appear in (0.0-1.0)

# Advanced options (usually not needed)
# script_dir: "/path/to/earlgrey/scripts"  # Auto-detected if installed via conda/mamba

# SLURM cluster settings (only used with --slurm flag)
slurm_partition: "normal.1000h"   # partition/queue to submit to (required with --slurm)
slurm_account: ""     # account string (leave empty if not required)
slurm_extra: ""       # any extra sbatch flags, e.g. "--constraint=avx2"
```

Submit the run, but check dry run first to confirm the DAG looks correct with 10 genomes and both modules enabled:

```bash
earlGreyParTEA \
  -c /data/toby/testDIR/0.1.6tests/8_full/manyGenomesTest.yaml \
  -t 32 \
  --slurm \
  --dry-run
```

This looks okay, let's run it for real!

```bash
earlGreyParTEA \
  -c /data/toby/testDIR/0.1.6tests/8_full/manyGenomesTest.yaml \
  -t 8 \
  --slurm
```

The many genomes test passed with no errors, and all expected outputs were produced. The shared/unique plots show the 10 species in the same order as the tree, and the BUSCO vs TE QC plot shows points coloured by dominant TE class. This confirms the pipeline is working correctly for v0.1.6.

---

## Release v0.1.7 Feature Updates

### Removal of `-norna` flag from all RepeatMasker calls

**Problem:** All RepeatMasker invocations in the pipeline carried the `-norna` flag, which suppresses masking of small RNA genes and related repeats (snRNA, scRNA, tRNA etc.). This was undesirable for a general-purpose TE annotation pipeline — RNA-derived repeats such as SINEs are legitimate TE families and excluding them from masking introduces a systematic gap in the annotation.

**Fix:** Removed `-norna` from every RepeatMasker call in the pipeline:

- `rules/lib_construct.smk` — removed from `repeatmasker` (species-database initial masking) and `repeatmasker_custom` (custom-library initial masking)
- `rules/annotate.smk` — removed from `repeatmasker_annotation` (final annotation masking against the curated TE library)

The warmup dummy run (`repeatmasker_warmup`) is not affected because it uses `-lib` with a short random sequence that will never match any annotated repeats regardless of `-norna`.

**Files changed:**
- `rules/lib_construct.smk`
- `rules/annotate.smk`

---

### Extended `repeatmasker_warmup` to pre-build species-specific library cache

**Background:** The v0.1.3 `repeatmasker_warmup` rule only pre-built the general RepeatMasker BLAST cache. When a taxon-specific library is requested (e.g. `repeatmasker_species: "7215"` for Drosophila), RepeatMasker must also build a species-specific BLAST database cache at `Libraries/CONS-Dfam_withRBRM_3.9/<species>/` the first time it runs. This cache build writes several files:

- `refineableHash.dat` — a Perl `Storable` hash mapping each TE family to its refinability status
- `speciesMeta.pm` — Perl module with metadata about the species library
- `specieslib` + BLAST DB files (`specieslib.nhr`, `.nin`, `.nsq`, etc.) — the BLAST-indexed library

If multiple genome jobs start simultaneously before this cache exists, each job tries to build it in parallel. All jobs create their own `<species>.working/` staging directory, but only one can win the final rename; all others fail immediately.

**Additional complication — incomplete cache from a previous OOM kill:** If a previous run was killed mid-build (e.g. by a SLURM cgroup OOM signal during `makeblastdb`), the `<species>/` directory may be left in a partially-built state containing the BLAST `.nhr`/`.nin`/`.nsq` files but lacking `refineableHash.dat` and `speciesMeta.pm`. RepeatMasker's cache-validation logic checks for `*.nhr` first; if found, it treats the cache as valid and tries to immediately retrieve `refineableHash.dat` with Perl `Storable::retrieve()`. This fails with a crash rather than triggering a rebuild, causing every parallel genome job to fail.

**Root cause of observed failures (April 2026 drosophila run):**
1. **April 24** — OOM kill during `makeblastdb` building the Drosophila species library (2,541 sequences). Left `7215/` with BLAST DB files only; no `refineableHash.dat`.
2. **April 30 re-run** — `7215/` had `.nhr` files → RepeatMasker's cache check passed → immediate crash on `retrieve("refineableHash.dat")` for all parallel jobs.

**Fix:** The `repeatmasker_warmup` rule now also pre-builds the species-specific cache before any parallel genome jobs start. The new logic:

1. Checks whether the species cache directory exists but is missing `refineableHash.dat` (incomplete from a previous aborted run). If so, **deletes the directory** to force a fresh build. RepeatMasker will not do this itself because it considers `.nhr` a sufficient validity indicator.
2. Removes any stale `<species>.working/` directory from a previous aborted parallel build.
3. Runs a single dummy RepeatMasker job (`-species <spec>` on a one-sequence FASTA) to trigger the full cache build — including `refineableHash.dat`, `speciesMeta.pm`, and the BLAST DB.
4. **Verifies** that `refineableHash.dat` now exists. If not, the warmup rule exits with a non-zero status and a clear error message, preventing the pipeline from proceeding with a broken cache.

The warmup `mem_mb` was increased from 4,000 to 32,000 MB and `runtime` from 30 to 120 minutes to accommodate the large species library build (the Drosophila `7215` library is ~2,500 sequences and requires a correspondingly large `makeblastdb` run).

**Files changed:**
- `rules/lib_construct.smk` — extended `repeatmasker_warmup` shell block; added `params: rep_spec=REPSPEC`; increased `mem_mb` and `runtime` resources

---

### RepeatModeler robustness for small genomes (`-genomeSampleSizeMax`)

**Problem:** RepeatModeler runs up to 6 rounds of repeat discovery. Round 1 uses RepeatScout and samples up to 40 Mbp; rounds 2–6 use RECON with progressively larger samples:

| Round | Tool | Sample size |
|-------|------|-------------|
| 1 | RepeatScout | 40 Mbp |
| 2 | RECON | 3 Mbp |
| 3 | RECON | 9 Mbp |
| 4 | RECON | 27 Mbp |
| 5 | RECON | 81 Mbp |
| 6 | RECON | 243 Mbp |

For genomes smaller than the round 6 default, RepeatModeler may attempt a round for which there is insufficient unmasked sequence (because earlier rounds have already masked much of the genome). When this happens it fails to write the round's `sampleDB-N.fa` file and crashes with:

```
FastaDB::compact - Error could not locate file .../round-N/sampleDB-N.fa!
 at .../RepeatModeler line 943.
```

This is a RepeatModeler bug — it should exit gracefully — but the pipeline must handle it.

**Complication — contigs < 40 kb are discarded:** RepeatModeler discards any contig shorter than 40 kb during sampling. Using the total genome size from the BLAST database (`blastdbcmd -info`) overestimates the sequence actually available. The correct value is the sum of all contig lengths ≥ 40 kb.

**Fix:** The `repeatmodeler` rule now:

1. Computes the **sampable genome size** — sum of contig lengths ≥ 40 kb — from the `.prep` FASTA using an awk one-liner (the `.prep` file is already an explicit rule input):

```bash
GENOME_SIZE=$(awk '/^>/{if(len>=40000)sum+=len; len=0; next}{len+=length($0)} \
    END{if(len>=40000)sum+=len; print sum+0}' {input.prep})
```

2. Selects the highest RECON round the genome can support by comparing the sampable size against **cumulative** sample totals across rounds (since a genome must have enough sequence for all rounds up to and including the cap, not just the final one):

| Sampable size | `-genomeSampleSizeMax` | Rounds run |
|---|---|---|
| ≥ 363 Mbp (3+9+27+81+243) | none (default) | 1–6 |
| ≥ 120 Mbp (3+9+27+81) | `81000000` | 1–5 |
| ≥ 39 Mbp (3+9+27) | `27000000` | 1–4 |
| ≥ 12 Mbp (3+9) | `9000000` | 1–3 |
| < 12 Mbp | `3000000` | 1–2 |

The `-genomeSampleSizeMax` value itself is set to the per-round size for the final allowed round (not the cumulative total), matching RepeatModeler's internal semantics for the flag.

**Example (Neuro73, ~38.7 Mbp sampable):** `39M > 38.7M`, so the cap is `9000000` → rounds 1–4 only. Without the fix RepeatModeler ran all 4 RECON rounds successfully then attempted round 5, could not build `sampleDB-4.fa` from the exhausted masked genome, and crashed.

**Files changed:**
- `rules/lib_construct.smk` — added `.prep` as an explicit input to `repeatmodeler`; awk sampable-size calculation; round-capping logic with cumulative thresholds

---

### Per-rule log files (all rules)

**Motivation:** Previously only a small number of rules had `log:` directives. All tool output (RepeatMasker, RepeatModeler, HELIANO, etc.) was printed to stdout/stderr and mixed with Snakemake's own progress messages. On long runs this made it hard to locate warnings or diagnose failures without scrolling through thousands of lines.

**Changes:** A `log:` directive was added to every rule across all six rule files:

| File | Rules modified |
|------|---------------|
| `rules/lib_construct.smk` | `repeatmasker_warmup`, `prep_genome`, `repeatmasker`, `repeatmasker_custom`, `extract_repeatmasker_library`, `build_db`, `repeatmodeler`, `testrainer` |
| `rules/annotate.smk` | `repeatmasker_annotation`, `heliano_detection`, `merge_repeats`, `generate_summary_charts`, `calculate_divergence`, `sweep_up_files`, `generate_softmasked_genome` |
| `rules/busco_phylo.smk` | `fetch_busco_db`, `run_busco`, `busco_summary_table`, `extract_busco_aa`, `align_busco_gene`, `create_supermatrix`, `run_fasttree`, `busco_completeness_phylo`, `busco_te_qc` |
| `rules/saturation.smk` | `saturation_plot` |
| `rules/clustering.smk` | `cluster_all_species` |
| `rules/shared_unique_content.smk` | `shared_unique_plot`, `shared_unique_plot_phylo`, `shared_unique_pa_plot`, `shared_unique_pa_plot_phylo` |

**Log file placement:** Log files sit alongside their output directories, using the same wildcard structure as the rule's `output:` block. Example paths:

- `{outdir}/{species}_EarlGrey/{species}_RepeatModeler/{species}.repeatmodeler.log`
- `{outdir}/{species}_EarlGrey/{species}_RepeatMasker/{species}.repeatmasker.log`
- `{outdir}/{species}_EarlGrey/{species}_heliano/{species}.heliano_detection.log`
- `{OUTDIR}/combinedLibraries/cluster_all_species.log`
- `{OUTDIR}/buscoPhylo/fasttree.log`

**Shell output redirection strategy:**

- `shell:` rules: `exec > {log} 2>&1` as the first line redirects all subsequent stdout and stderr (including subprocesses) to the log file.
- `script:` rules: Snakemake automatically redirects stderr to the `log:` file; no extra line is needed.
- `run:` rules: each `shell()` call appends `>> str(log) + " 2>&1"`.
- **Special case — `run_fasttree`:** FastTree writes the newick tree to stdout and diagnostics to stderr, so `FastTree ... > {output.tree} 2>> {log}` is used instead of `exec > {log} 2>&1`.

**Effect:** After this change only Snakemake's job submission/completion messages and explicit `print()` calls in `on_start_functions.py` appear on the terminal. All tool output is captured in per-rule log files.

**Files changed:**
- `rules/lib_construct.smk`, `rules/annotate.smk`, `rules/busco_phylo.smk`, `rules/saturation.smk`, `rules/clustering.smk`, `rules/shared_unique_content.smk` — `log:` directives and redirection added to all rules

---

### Bug fix: wildcard mismatch in `clustering.smk`

**Error observed when the log directive was first added:**

```
WorkflowError: Not all output, log and benchmark files of rule cluster_all_species
contain the same wildcards. This is a Snakemake requirement.
```

**Root cause:** The `log:` directive was initially written as a Python f-string:

```python
log: f"{OUTDIR}/combinedLibraries/cluster_all_species.log"
```

Because f-strings are evaluated at Python parse time, this produces a literal string that Snakemake's wildcard system cannot see. The rule's `output:` uses the `{outdir}` Snakemake wildcard, so Snakemake's validator correctly rejects the mismatch.

**Fix:** Changed to a plain wildcard string:

```python
log: "{outdir}/combinedLibraries/cluster_all_species.log"
```

**General principle:** `log:` must use the same wildcard form as `output:`. If `output:` uses Snakemake wildcards (`{outdir}`, `{species}`, etc.), so must `log:`. Python f-strings in rule fields are only appropriate when `output:` is also an f-string (path fully determined at parse time from global variables like `OUTDIR`).

**Files changed:**
- `rules/clustering.smk` — `log:` changed from f-string to Snakemake wildcard string

---

### Verification checklist for v0.1.7

- [ ] Run pipeline with `repeatmasker_species` set on a species with a large library (e.g. Drosophila `7215`). Confirm warmup builds a complete `7215/` cache including `refineableHash.dat` before any parallel genome jobs start.
- [ ] Manually corrupt the `7215/` cache (delete `refineableHash.dat`, leave `.nhr` files). Re-run. Confirm warmup detects the incomplete cache, removes the directory, and rebuilds it cleanly.
- [ ] Run pipeline with a genome < 81 Mbp. Confirm RepeatModeler receives `-genomeSampleSizeMax <size>` in the job log.
- [ ] Run pipeline with a genome > 81 Mbp. Confirm RepeatModeler receives no extra flag.
- [ ] Confirm no RNA-family annotations appear as masked in RepeatMasker output that would previously have been suppressed by `-norna` (e.g. check for `RNA` or `srpRNA` class entries in `.out` files).
- [ ] After a completed run, confirm per-rule log files exist alongside each output directory (e.g. `{species}.repeatmodeler.log`, `{species}.repeatmasker.log`, `{species}.heliano_detection.log`).
- [ ] Confirm terminal output during a run shows only Snakemake progress lines and startup messages — no RepeatMasker/RepeatModeler/HELIANO stdout.
- [ ] Confirm `{OUTDIR}/combinedLibraries/cluster_all_species.log` is created on a multi-genome run (validates the wildcard-string fix).
- [ ] Local mode regression: confirm all changes do not break a standard local run on a previously working dataset.

## Release v0.1.8 Feature Updates

### CD-HIT length-difference cutoff for clustering (`clustering_length_diff`)

**Problem:** The `cd-hit-est` command was invoked with `-aS` (alignment coverage of the shorter sequence) and `-G 0` (local alignment mode) but no constraint on the *length ratio* between sequences. Under these settings a small partial TE consensus (e.g. 1 000 nt) could be placed in the same cluster as a full-length copy (e.g. 20 000 nt) provided ~80 % of the short sequence aligned somewhere within the long one. Clusters spanning a 10–20× length range can merge biologically distinct elements (e.g. an LTR solo with a complete LTR retrotransposon) and produce misleading cluster representatives.

**Fix:** A new `clustering_length_diff` parameter is exposed in `config.yaml` (default `0.5`) and passed to `cd-hit-est` as the `-s` flag. This requires the shorter sequence to be at least this fraction of the longer sequence's length before clustering is allowed. At the default value of `0.5`, sequences must be within a 2× length ratio — removing the most extreme mismatches while still permitting partial elements from the same family to cluster.

Set to `0.0` to restore the previous behaviour (no length restriction).

**Files changed:**
- `config/config.yaml` — new `clustering_length_diff: 0.5` parameter
- `rules/clustering.smk` — `cluster_length_diff` param added; `-s {params.cluster_length_diff}` added to `cd-hit-est` call
- `scripts/on_start_functions.py` — default added to `validate_parameters`; startup message updated to include `length_diff`
- `earlGreyParTEA`, `earlGreyParTEA_LibConstruct`, `earlGreyParTEA_AnnotationOnly` — `clustering_length_diff: 0.5` added to generated config templates
- `README.md` — v0.1.8 section updated; Clustering Options docs updated

---

### CD-HIT long-sequence coverage filter (`clustering_coverage_long` / `-aL`)

**Problem:** `cd-hit-est` was run with `-G 0 -aS 0.8` (80% of the *shorter* sequence must align). This allowed a short sequence (e.g. 1 092 nt) to be absorbed into a cluster whose representative was an order of magnitude longer (e.g. 20 383 nt), provided ~874 nt aligned somewhere in the longer sequence. Inspection of the resulting `.clstr` files revealed that some very long representatives appear to be chimeric or over-extended consensus sequences — they contain sub-regions with similarity to multiple distinct short elements, causing those short elements to be incorrectly merged into a single cluster.

**Analysis:** Examining alignment coordinates in the `.clstr` format reveals the mechanism. In a chimeric representative of 20 383 nt, different members align to completely non-overlapping windows (e.g. positions 1–13 278 and 16 170–17 256), confirming the representative spans at least two distinct elements. The `-aS` parameter cannot guard against this because a short sequence's coverage of itself is always high regardless of what fraction of the long representative it covers.

**Fix:** The `-aL` flag requires the alignment to cover at least a given fraction of the *longer* sequence. Because the alignment cannot be longer than the shorter sequence, this imposes an implicit length-ratio constraint:

$$\text{alignment} \geq \text{aL} \times \text{longer} \quad\Rightarrow\quad \frac{\text{shorter}}{\text{longer}} \geq \text{aL}$$

| `clustering_coverage_long` | Max tolerated length ratio | Example effect |
|---|---|---|
| 0.0 (default, disabled) | no limit | previous behaviour |
| 0.5 | ~2× | 7 779 nt cannot cluster with 20 383 nt representative |
| 0.75 | ~1.33× | both 7 779 nt and 14 766 nt excluded from 20 383 nt cluster |
| 0.8 | ~1.25× | strict; sequences must be nearly the same length |

Recommended value: `0.75`. Default is `0.0` for backward compatibility.

**Files changed:**
- `config/config.yaml` — new `clustering_coverage_long: 0.0` parameter
- `rules/clustering.smk` — `cluster_coverage_long` param added; `-aL {params.cluster_coverage_long}` added to `cd-hit-est` call
- `scripts/on_start_functions.py` — default added; startup message extended to report `aL` status
- `earlGreyParTEA`, `earlGreyParTEA_LibConstruct`, `earlGreyParTEA_AnnotationOnly` — `clustering_coverage_long: 0.0` added to generated config templates
- `README.md` — Clustering Options docs updated

---

### Dynamic discovery of RepeatMasker species-library cache directory

**Problem:** The `repeatmasker_warmup` rule hardcoded the species-library BLAST cache parent directory as `$RM_SHARE/Libraries/CONS-Dfam_withRBRM_3.9`. This path is only present when the RepeatMasker installation includes both Dfam **and** RepBase RepeatMasker edition (RBRMSK). Users who configured RepeatMasker with Dfam only (or a future Dfam version with a different suffix) have a differently-named directory (e.g. `CONS-Dfam_3.9`). In those environments the warmup checked and rebuilt the cache in the wrong location — the real cache directory was never validated before the parallel genome jobs started.

**Fix:** The hardcoded `CACHE_PARENT` assignment is replaced by a `find` call that discovers the actual `CONS-*` directory at runtime:

```bash
CACHE_PARENT=$(find "$RM_SHARE/Libraries" -maxdepth 1 -type d -name "CONS-*" 2>/dev/null | head -n 1)
```

This executes after the general-library warmup (which ensures the `Libraries/` directory is fully populated), so the `CONS-*` directory will exist by the time `find` runs if any species library is present. If no matching directory is found the warmup prints a clear warning and exits gracefully (exit 0) rather than attempting to validate a cache inside a non-existent parent.

**Files changed:**
- `rules/lib_construct.smk` — hardcoded `CACHE_PARENT` path replaced with dynamic `find` discovery; added graceful exit with warning if no `CONS-*` directory is found

---

### Post-clustering chimera detection and cluster splitting (`split_chimeras`)

**Background:** Even after applying `-s` (length-difference) and `-aL` (long-sequence coverage) constraints at clustering time, some chimeric cluster representatives can survive. A chimeric representative is a consensus sequence built from two or more distinct TE families joined end-to-end during EarlGrey's iterative BEAT process. When such a representative acts as a cluster hub, biologically unrelated short elements map to opposite ends of it and end up incorrectly merged into a single cluster.

**Detection algorithm:**
For each cluster with ≥ `chimera_min_members` non-representative members:
1. Extract `(rep_start, rep_end)` alignment coordinates for each member from the `.clstr` file (available because `cd-hit-est` is run with `-d 0`).
2. Build an overlap graph: two members share an edge if their alignment windows on the representative overlap by ≥ `chimera_overlap_min` nucleotides (default 50 nt).
3. Find connected components via BFS. If ≥ 2 components each spanning ≥ `chimera_min_component_span` fraction of the representative's length are found, the cluster is flagged as chimeric.
4. Assign a **chimera score** = largest inter-component gap (nt) / representative length. Higher scores indicate a clearer structural break.

**Splitting:**
For each confirmed chimeric cluster:
- The original representative is written to the output FASTA with a `_CHIMERA` suffix appended to its base name, retaining its `#classification` tag for traceability.
- For each component (sorted by leftmost alignment position on the representative), the **longest member** is selected as the new cluster representative. Its sequence and **full original FASTA header** are taken directly from the pre-clustering combined FASTA (`combined_all_species.fa`), so all species prefixes and classification tags are preserved exactly as they were entered into clustering.

**Snakemake integration:**
- `cluster_all_species` now declares `combined_all_species.fa` as a named `temp()` output rather than deleting it with `rm -f`. Snakemake manages its deletion automatically once all consuming rules finish.
- The `split_chimeras` rule is only defined (and only appears in the DAG) when `split_chimeras: true` AND `skip_clustering: false`.
- `rules/annotate.smk` resolves the annotation library at parse time: `combined_all_species.chimera_split.fa` when `split_chimeras: true`, otherwise `combined_all_species.clstrd.fa`. The original clustered FASTA is always preserved.
- `combined_all_species.chimera_split.fa` and `chimera_detection_summary.tsv` are requested by `rule all` for `full` and `libconstruct` pipeline modes when the feature is active.

**New config parameters (all backward-compatible defaults):**
```yaml
split_chimeras: false           # enable post-clustering chimera detection and splitting
chimera_overlap_min: 50         # min nt overlap between member alignment windows to be
                                 #   in the same component (lower = more sensitive)
chimera_min_members: 3          # min non-representative members to test a cluster
chimera_min_component_span: 0.1 # each component must span >= this fraction of rep length
```

**Output files** (written to `{output_dir}/combinedLibraries/`):
- `combined_all_species.chimera_split.fa` — modified library used for annotation; chimeric reps relabelled `_CHIMERA`, replaced by component representatives with original headers
- `chimera_detection_summary.tsv` — per-cluster table with columns: `cluster_idx`, `representative`, `rep_length`, `n_members`, `is_chimeric`, `n_components`, `component_sizes`, `chimera_score`, `component_rep_names`

**Files changed:**
- `scripts/split_chimeras.py` — new script; detection, overlap graph, BFS, and FASTA output; guarded with `if 'snakemake' in dir(): main()` so the module is importable for testing
- `rules/clustering.smk` — `combined_fa=temp(...)` added to `cluster_all_species` outputs; `rm -f` lines removed; `split_chimeras` rule appended (conditionally defined)
- `rules/annotate.smk` — `ANNOTATION_LIBRARY` variable resolves to `chimera_split.fa` or `clstrd.fa` based on config
- `Snakefile` — chimera outputs added to `rule all` for `libconstruct` and `full` modes
- `config/config.yaml`, `earlGreyParTEA`, `earlGreyParTEA_LibConstruct`, `earlGreyParTEA_AnnotationOnly` — four new chimera params added to config templates
- `scripts/on_start_functions.py` — four new defaults added; startup message reports chimera detection status and params
- `README.md` — feature bullet added; new v0.1.8 subsection; new `### Chimera Detection Options` config section

---

### Verification checklist for v0.1.8

> **Key:** `[x]` = verified by automated tests (`tests/test_v018_features.py`);
> `[ ]` = still requires a real conda environment or pipeline run (see commands below).

- [ ] Run pipeline with `repeatmasker_species` set on an installation configured with Dfam **only** (no RepBase). Confirm warmup locates the `CONS-Dfam_*` directory and validates/builds the species cache correctly.
- [ ] Run pipeline on a standard installation with both Dfam and RepBase. Confirm existing behaviour is unchanged.
- [x] Confirm that if the `Libraries/` directory contains no `CONS-*` subdirectory the warmup prints the expected warning and the pipeline does not proceed with species-cache validation. *(TestCacheDiscovery::test_no_cons_dir_exits_zero_with_warning)*
- [x] Run clustering on a multi-genome dataset and confirm that sequences with a length ratio > 2× (e.g. ~1 000 nt vs ~20 000 nt) are no longer grouped in the same cluster. *(TestCdhitIntegration::test_length_diff_prevents_extreme_size_mismatch)*
- [x] Confirm that setting `clustering_length_diff: 0.0` in the config restores the previous behaviour (no length restriction), and that sequences of very different sizes can still be clustered. *(TestCdhitIntegration::test_length_diff_zero_allows_clustering)*
- [x] Confirm the startup summary prints `length_diff:` alongside `identity:` and `coverage:` when clustering is enabled. *(TestStartupMessages::test_length_diff_reported_in_message)*
- [x] Confirm `--generate-config` output from all three wrapper scripts includes `clustering_length_diff: 0.5`. *(TestGenerateConfig — full + libconstruct wrappers)*
- [x] Set `clustering_coverage_long: 0.75` and re-run clustering on a dataset with known chimeric clusters (e.g. arabidopsis combinedLibraries). Confirm sequences aligning to <75% of the representative's length are no longer placed in the same cluster. *(TestCdhitIntegration::test_aL_prevents_short_in_long_cluster)*
- [x] Confirm `clustering_coverage_long: 0.0` (default) restores the original `-aL 0.0` behaviour. *(TestCdhitIntegration::test_aL_zero_does_not_restrict)*
- [x] Confirm the startup message correctly reports `aL: disabled` when `clustering_coverage_long: 0.0` and reports the value when non-zero. *(TestStartupMessages::test_aL_disabled_message_when_zero + test_aL_value_reported_when_nonzero)*
- [x] Confirm `--generate-config` output from all three wrapper scripts includes `clustering_coverage_long: 0.0`. *(TestGenerateConfig — full + libconstruct wrappers)*
- [ ] Enable `split_chimeras: true` on a multi-genome dataset and confirm `combined_all_species.chimera_split.fa` and `chimera_detection_summary.tsv` are produced.
- [x] Verify that chimeric representatives appear with `_CHIMERA` suffix in `chimera_split.fa` and that component representatives carry their original FASTA headers unchanged. *(TestMainIntegration::test_chimeric_rep_labelled + test_component_reps_have_original_headers)*
- [ ] Confirm that when `split_chimeras: true`, downstream RepeatMasker annotation uses `chimera_split.fa` (check Snakemake DAG or log for the correct input path).
- [ ] Confirm `split_chimeras: false` (default) leaves pipeline behaviour completely unchanged and `clstrd.fa` is used for annotation.
- [ ] Confirm that `split_chimeras: true` with `skip_clustering: true` produces a clear error or is silently ignored (rule is not defined in DAG).
- [x] Run the unit-importable test (`python3 -c "from scripts.split_chimeras import parse_clstr, detect_chimera"`) and confirm no `NameError` is raised. *(TestModuleImportability)*

---

### Commands for remaining manual verifications

#### Items 4, 5, 8, 9 — cd-hit-est behavioural tests (run inside conda env)

```bash
# Requires cd-hit-est on PATH (automatically available in the conda env)
pytest tests/test_v018_features.py -v -m integration
```

Expected: 4 tests pass (`test_length_diff_prevents_extreme_size_mismatch`, `test_length_diff_zero_allows_clustering`, `test_aL_prevents_short_in_long_cluster`, `test_aL_zero_does_not_restrict`).

These passed!

#### Items 1, 2 — RepeatMasker warmup with a real species library

Items 1 and 2 require a full pipeline run on a machine with RepeatMasker and a real Dfam/RepBase species library installed. Add `repeatmasker_species` to config and run:

```bash
cd /data/toby/testDIR/drosophila_testSet

earlGreyParTEA --generate-config test12_018.yaml --genome-dir /data/toby/testDIR/drosophila_testSet --output-dir /data/toby/testDIR/

earlGreyParTEA -c /data/toby/testDIR/drosophila_testSet/test12_018.yaml --threads 16

# After the warmup step completes, confirm the cache was built:
RM_SHARE=$(which RepeatMasker | sed 's|/bin/RepeatMasker$|/share/RepeatMasker|')
find "$RM_SHARE/Libraries" -maxdepth 1 -type d -name "CONS-*"
# /data/toby/miniforge3/envs/partea_018/share/RepeatMasker/Libraries/CONS-Dfam_withRBRM_3.9

find $RM_SHARE/Libraries/CONS-*/7215 -name "refineableHash.dat"
# /data/toby/miniforge3/envs/partea_018/share/RepeatMasker/Libraries/CONS-Dfam_withRBRM_3.9/7215/refineableHash.dat
```

For item 2 (standard installation, unchanged behaviour) simply run the same config and confirm the pipeline completes.normally.

This run worked perfectly! The warmup built the `7215/` cache with `refineableHash.dat` as expected, and the pipeline completed without issue.

#### Item 12 — Full `split_chimeras` pipeline run

```bash
# Add to config.yaml:
#   split_chimeras: true
#   chimera_overlap_min: 50
#   chimera_min_members: 3
#   chimera_min_component_span: 0.1

earlGreyParTEA --generate-config /data/toby/testDIR/test_splitChimeras_018/test_splitChimeras_018.yaml --genome-dir /data/toby/testDIR/drosophila_testSet --output-dir /data/toby/testDIR/test_splitChimeras_018

earlGreyParTEA -c /data/toby/testDIR/test_splitChimeras_018/test_splitChimeras_018.yaml --threads 16

# Check both output files were created:
ls <output_dir>/combinedLibraries/combined_all_species.chimera_split.fa
ls <output_dir>/combinedLibraries/chimera_detection_summary.tsv

# Quick sanity check — count chimeric entries:
grep -c "_CHIMERA" <output_dir>/combinedLibraries/combined_all_species.chimera_split.fa
```

No chimeras were detected in this run. This is expected: see **"When chimeric TEs are actually found"** below.

---

#### When chimeric TEs are actually found

**Why the drosophila test produced no chimeras**

Two factors suppress chimera detection in this run:

1. **`clustering_coverage_long` (the `-aL` flag) actively prevents chimeric clustering from forming.** When `-aL > 0`, a short sequence can only cluster with a much longer representative if the alignment covers ≥ `aL` fraction of the representative's length. Members from a second TE family that align only to a narrow window of a long chimeric representative (e.g. covering only 10–30% of the rep) do not meet this threshold and are not placed in the cluster at all. With those members absent, the two populations of alignment windows needed to trigger chimera detection never appear in the `.clstr` file.

2. **The drosophila test set is a small, well-annotated dataset.** Chimeric BEAT consensus sequences most commonly arise from larger, less-characterised genomes where EarlGrey builds consensus sequences across complex TE landscapes.

**The paradox: tighter clustering filters → fewer chimeras detectable**

| Parameter | Effect on chimera detection |
|---|---|
| `clustering_coverage_long: 0.0` (default, disabled) | Most permissive — sequences of any length can cluster, including short members that reveal chimeric structure |
| `clustering_coverage_long: 0.75` | Short sequences aligning to < 75 % of the representative are excluded → fewer members per chimeric cluster → harder to detect |
| `clustering_length_diff: 0.0` (disabled) | No length ratio constraint → extreme-length pairs can cluster |
| `clustering_length_diff: 0.5` (default) | Sequences must be within 2× length ratio |

The `-aL` and `-s` filters primarily *prevent* chimeric clustering from occurring; `split_chimeras` catches chimeras that survive despite those filters. Running clustering with `-aL 0.0` (disabled) will allow more members to enter each cluster and makes chimeras more detectable, at the cost of allowing more extreme length mismatches.

**Conditions most likely to produce detectable chimeras**

- `clustering_coverage_long: 0.0` (disabled)
- Large, complex genomes (e.g. plants, polyploids, large arthropods)
- Datasets with many LTR retrotransposons (LTR termini are shared across families and seed chimeric BEAT consensus building)
- Many genomes pooled: more opportunities for a chimeric representative from one genome to absorb members from others

**Standalone chimera detection test (fixture files)**

Test fixture files in `tests/fixtures/` contain a synthetic chimeric cluster:

| File | Contents |
|---|---|
| `chimera_test.clstr` | 2 clusters; cluster 0 has members aligning to pos 1–500 (DNA/hAT) and 1400–2000 (LINE/L1) on a 2000 nt rep |
| `chimera_test_clustered.fa` | Cluster representative sequences |
| `chimera_test_combined.fa` | All pre-clustering sequences with original headers |

Expected output of `split_chimeras.py` on these fixtures:
- `rep_chimeric_CHIMERA#Unknown` — original chimeric representative labelled and retained
- `>sp1_hAT#DNA/hAT EarlGrey_annotation` — component 1 rep (longest DNA/hAT member), original header preserved
- `>sp1_L1#LINE/L1 EarlGrey_annotation` — component 2 rep (longest LINE/L1 member), original header preserved
- `sp1_Tc1#DNA/TcMar-Tc1` — clean non-chimeric cluster passes through unchanged
- Chimera score = 0.45 (900 nt gap / 2000 nt rep length)

Run the fixture test via the existing unit test suite:

```bash
cd /data/toby/EarlGreyParTEA
conda run -n partea_018 pytest tests/test_v018_features.py::TestMainIntegration -v
```

Or run the script directly against the fixture files:

```bash
cd /data/toby/EarlGreyParTEA
conda run -n partea_018 python3 - <<'PYEOF'
import sys; sys.path.insert(0, 'scripts')
import split_chimeras as sc
from types import SimpleNamespace

FIXTURES = 'tests/fixtures'
sc.__dict__['snakemake'] = SimpleNamespace(
    input=SimpleNamespace(
        clstr=f'{FIXTURES}/chimera_test.clstr',
        clustered_fa=f'{FIXTURES}/chimera_test_clustered.fa',
        combined_fa=f'{FIXTURES}/chimera_test_combined.fa',
    ),
    output=SimpleNamespace(
        fasta='/tmp/chimera_split.fa',
        summary='/tmp/chimera_detection_summary.tsv',
    ),
    log=['/tmp/chimera_test.log'],
    params=SimpleNamespace(overlap_min=50, min_members=3, min_component_span=0.1),
)
sc.main()
for line in open('/tmp/chimera_split.fa'):
    if line.startswith('>'): print(line.rstrip())
print()
print(open('/tmp/chimera_detection_summary.tsv').read())
PYEOF
```

---

#### Item 14 — Confirm annotation uses `chimera_split.fa` when `split_chimeras: true`

```bash
# Dry-run with split_chimeras: true — inspect the planned annotation job input:
earlGreyParTEA -c <your_config_split_chimeras_true.yaml> --dry-run 2>&1 \
  | grep -A5 "repeatmasker_annotation"
# → input line should reference combined_all_species.chimera_split.fa

# Or after a real run, check the annotation log:
grep "chimera_split\|clstrd" <output_dir>/<species>_EarlGrey/<species>_RepeatMasker/<species>.repeatmasker_annotation.log | head -5
```

---

#### Item 15 — Confirm `split_chimeras: false` leaves pipeline unchanged (`clstrd.fa` used)

```bash
# Dry-run with split_chimeras: false (the default):
earlGreyParTEA -c <your_config_default.yaml> --dry-run 2>&1 \
  | grep -A5 "repeatmasker_annotation"
# → input line should reference combined_all_species.clstrd.fa

# Confirm split_chimeras rule does NOT appear in the DAG:
earlGreyParTEA -c <your_config_default.yaml> --dry-run 2>&1 \
  | grep "split_chimeras"
# → no output expected
```

---

#### Item 16 — Confirm `split_chimeras: true` + `skip_clustering: true` → rule absent from DAG

```bash
# Add to your config.yaml:
#   split_chimeras: true
#   skip_clustering: true

earlGreyParTEA -c <your_config_skip_clust.yaml> --dry-run 2>&1 \
  | grep "split_chimeras"
# → no output expected (rule is conditionally defined only when skip_clustering: false)
```

## Release v0.1.9 Feature Updates

### LSF cluster submission via `--lsf` (Experimental)

**Background:** A user request to support LSF (IBM Spectrum LSF / Platform LSF) cluster environments in addition to SLURM. LSF users at HPC centres cannot use the `--slurm` mode, leaving them with only the local execution path. The implementation follows the same pattern as the existing SLURM integration: each Snakemake rule is submitted as an individual cluster job, with CPUs, memory, and wall-time derived from the rule's existing `resources:` block.

---

#### Approach: community LSF executor plugin

The implementation uses `snakemake-executor-plugin-lsf`, activated via `--executor lsf`. This is a community-maintained plugin (author: Brian Fulton-Howard) available from the Snakemake plugin catalog. It is not part of the official Snakemake organisation, which is why the feature is flagged as **experimental** throughout.

Resource mapping from rule `resources:` block to `bsub` flags:
- `threads` → `-n` (number of cores per job)
- `mem_mb` → `-R rusage[mem=<mem_mb/threads>]` (per-core memory by default; set `SNAKEMAKE_LSF_MEMFMT=perjob` for per-job total)
- `runtime` → `-W <minutes>` (wall-time limit)
- `lsf_queue` → `-q` (queue/partition)
- `lsf_project` → `-P` (project string)
- `lsf_extra` → appended verbatim to the `bsub` command

No changes to any rule's `resources:` block are required — the same declarations used for SLURM are reused directly.

---

#### New CLI flags (all three entry points)

| Flag | Description |
|------|-------------|
| `--lsf` | Enable LSF submission mode |
| `--lsf-queue QUEUE` | LSF queue (required; can also be set as `lsf_queue:` in config) |
| `--lsf-jobs N` | Max concurrent LSF jobs (default: genomes × 3, capped at 200) |
| `--lsf-project PROJ` | LSF project string (optional) |
| `--lsf-extra "FLAGS"` | Extra raw `bsub` flags (e.g. `"-R 'select[type==X86_64]'"`) |

`--slurm` and `--lsf` are mutually exclusive. Using both simultaneously prints an error and exits.

---

#### Snakemake invocation in LSF mode

```bash
DEFAULT_RESOURCES=("lsf_queue=$LSF_QUEUE")
[ -n "$LSF_PROJECT" ] && DEFAULT_RESOURCES+=("lsf_project=$LSF_PROJECT")
[ -n "$LSF_EXTRA" ]   && DEFAULT_RESOURCES+=("lsf_extra=$LSF_EXTRA")

snakemake \
    --snakefile "$SNAKEFILE" \
    --configfile "$CONFIG" \
    --config lsf_mode=true \
    --executor lsf \
    --default-resources "${DEFAULT_RESOURCES[@]}" \
    --jobs "$LSF_JOBS" \
    --cores "$THREADS" \
    --latency-wait 60 \
    --retries 1 \
    $DRY_RUN $UNLOCK $RERUN
```

Key differences from SLURM invocation:
- `--executor lsf` instead of `--executor slurm`
- `--config lsf_mode=true` instead of `slurm_mode=true`
- Default resources use `lsf_queue`, `lsf_project`, `lsf_extra` instead of `slurm_partition`, `slurm_account`, `slurm_extra`

---

#### Thread lambda modification

The per-genome thread lambdas in `rules/lib_construct.smk`, `rules/annotate.smk`, and `rules/busco_phylo.smk` were already updated for SLURM mode with `config.get("slurm_mode", False)`. The condition is extended to also check `lsf_mode`:

```python
# Before
if config.get("slurm_mode", False)

# After
if config.get("slurm_mode", False) or config.get("lsf_mode", False)
```

This ensures each job gets the full `-t` thread count in LSF mode rather than the divided-across-genomes value used in local mode.

**Rules modified:** `repeatmasker`, `repeatmasker_custom`, `repeatmodeler`, `testrainer` (in `lib_construct.smk`); `repeatmasker_annotation` ×2, `heliano_detection`, merge-area rules (in `annotate.smk`); `run_busco` (in `busco_phylo.smk`).

---

#### Memory note: per-core vs per-job

`snakemake-executor-plugin-lsf` divides `mem_mb` by `threads` to produce a per-core memory request by default. Most LSF clusters expect per-core values in `-R rusage[mem=...]`. Users whose cluster accepts per-job totals can override with:

```bash
export SNAKEMAKE_LSF_MEMFMT=perjob
```

This is documented in the README and printed as part of the experimental warning at runtime.

---

#### Config additions

Three new optional keys added to `config/config.yaml` and all three `generate_config.py` templates. They are only used when `--lsf` is passed:

```yaml
# LSF cluster settings (only used with --lsf flag) [EXPERIMENTAL]
lsf_queue: ""         # queue to submit to (required when using --lsf; can be set here instead of --lsf-queue)
lsf_project: ""       # project string (leave empty if not required)
lsf_extra: ""         # any extra bsub flags, e.g. "-R 'select[type==X86_64]'"
```

Three silent defaults added to `validate_parameters()` in `scripts/on_start_functions.py`:
```python
'lsf_queue':   ("", None),
'lsf_project': ("", None),
'lsf_extra':   ("", None),
```

---

#### Inline heredoc fix in AnnotationOnly and LibConstruct wrappers

During implementation it was discovered that the inline `generate_example_config()` heredocs in `earlGreyParTEA_AnnotationOnly` and `earlGreyParTEA_LibConstruct` were missing the SLURM block entirely (the SLURM block had been added to the main `earlGreyParTEA` wrapper but not propagated to the other two). Both wrappers now include both a SLURM block and an LSF block in their heredoc templates, matching the main wrapper and the `generate_config.py` templates.

---

#### Files changed

- `conda/meta.yaml` — added `snakemake-executor-plugin-lsf` as a run dependency alongside `snakemake-executor-plugin-slurm`
- `earlGreyParTEA` — added `--lsf`, `--lsf-queue`, `--lsf-jobs`, `--lsf-project`, `--lsf-extra` flags; mutual exclusion guard; `elif $LSF` build block; LSF block added to inline heredoc
- `earlGreyParTEA_AnnotationOnly` — same LSF changes + SLURM block added to inline heredoc (previously missing)
- `earlGreyParTEA_LibConstruct` — same LSF changes + SLURM block added to inline heredoc (previously missing)
- `config/config.yaml` — LSF block added
- `scripts/generate_config.py` — LSF block added to `FULL_TEMPLATE`, `LIBCONSTRUCT_TEMPLATE`, `ANNOTATE_TEMPLATE`
- `scripts/on_start_functions.py` — three LSF defaults added to `validate_parameters()`
- `rules/lib_construct.smk` — thread lambdas updated to check `lsf_mode` (4 rules)
- `rules/annotate.smk` — thread lambdas updated to check `lsf_mode` (4 rules)
- `rules/busco_phylo.smk` — thread lambda updated to check `lsf_mode` (1 rule)
- `README.md` — new `## 🖥️ LSF Cluster Submission (Experimental)` section; LSF feature bullet added (chronological); `## Requirements` updated with LSF plugin entry; LSF entry added to TOC

---

### Verification checklist for v0.1.9

> All items require an LSF cluster to fully verify. Structural and unit-level checks can be done locally.

- [ ] `--lsf --lsf-queue myqueue` produces the expected Snakemake invocation with `--executor lsf` (check via `--dry-run`)
- [ ] `--slurm` and `--lsf` together print error and exit with non-zero status
- [ ] `--generate-config` output from all three wrappers includes both `slurm_partition`, `slurm_account`, `slurm_extra` AND `lsf_queue`, `lsf_project`, `lsf_extra` blocks
- [ ] `validate_parameters()` does not raise on a config with no LSF keys (silent defaults applied)
- [ ] On a real LSF cluster: submit a small dry run and confirm `bsub` jobs appear in `bjobs` with correct queue, CPUs, and memory
- [ ] Retry logic: kill a job manually and confirm Snakemake resubmits at 2× memory
- [ ] `SNAKEMAKE_LSF_MEMFMT=perjob` changes the memory request format as documented

### Commands for local pre-release checks

```bash
# Confirm mutual exclusion guard triggers correctly
./earlGreyParTEA --slurm --lsf -c config/config.yaml 2>&1 | grep -i "cannot\|error\|exclusive"

# Dry run with --lsf to inspect Snakemake invocation
./earlGreyParTEA -c /data/toby/testDIR/test.yaml -t 4 --lsf --lsf-queue testqueue --dry-run 2>&1 | head -40

# Confirm generate-config includes LSF block
./earlGreyParTEA --generate-config /tmp/test_v019.yaml
grep -A5 "lsf_queue" /tmp/test_v019.yaml

# Confirm AnnotationOnly and LibConstruct generate-config also include both SLURM and LSF blocks
./earlGreyParTEA_AnnotationOnly --generate-config /tmp/test_v019_anno.yaml
grep "slurm_partition\|lsf_queue" /tmp/test_v019_anno.yaml

./earlGreyParTEA_LibConstruct --generate-config /tmp/test_v019_lib.yaml
grep "slurm_partition\|lsf_queue" /tmp/test_v019_lib.yaml
```

### Bump version to 0.1.9

- [x] `earlGreyParTEA` — `VERSION="0.1.9"`
- [x] `earlGreyParTEA_LibConstruct` — `VERSION="0.1.9"`
- [x] `earlGreyParTEA_AnnotationOnly` — `VERSION="0.1.9"`
- [x] `conda/meta.yaml` — `{% set version = "0.1.9" %}`
- [x] `.github/workflows/conda-release.yml` — `default: '0.1.9'`
- [x] `README.md` — version bump in TOC link; new `## 🆕 Changes in Latest Release (v0.1.9)` section

This is now working as expected. I have a bug with the plots where the bars do not line up with the phylogenetic tree. This needs to be resolved before release. 


---

## Release v0.2.0 — Dfam 4.0 / FamDB 3.0.0 / RepeatMasker 4.2.4 compatibility

### Background

EarlGrey 7.3.0 upgraded its toolchain to RepeatMasker 4.2.4, RMBlast 2.17.0, RepeatModeler
2.0.9, and FamDB 3.0.0 (Dfam 4.0). The most significant structural change is that FamDB is
now a standalone conda package (`share/famdb-3.0.0/`) rather than embedded inside RepeatMasker's
`share/RepeatMasker/Libraries/famdb/` directory. RepeatMasker itself is pre-configured by the
conda post-install hook and points to the standalone famdb via `FAMDB_DIR` in
`RepeatMaskerConfig.pm`. No `perl ./configure` step is needed by the user.

The Dfam 4.0 HDF5 partition files are now named `dfam40.*.h5` (previously `dfam39_full.*.h5`).
The species-library BLAST cache directory is now `CONS-Dfam_4.0/` (previously `CONS-Dfam_3.9/`
or `CONS-Dfam_withRBRM_3.9/`), but the existing dynamic `CONS-*` discovery in the
`repeatmasker_warmup` rule already handles this without any change.

The new user-facing setup flow is simply: install via conda/mamba, then run
`download_dfam.py` (an interactive download tool that ships with the `famdb` package) to
fetch the desired Dfam 4.0 partitions.

### Root causes of breakage

1. **`extract_repeatmasker_library` shell block** — derived famdb path as
   `share/RepeatMasker/Libraries/famdb/`, which in the new environment contains only
   `RMRBMeta.embl` / `RMRB_DUP.txt` (not the HDF5 files). The `famdb.py -i $libpath`
   call would fail or return an empty library.

2. **`validate_parameters` EarlGrey configuration check** — also derived the old path,
   then checked for `dfam39_full.*.h5` files and the `.earlgrey.config.complete` marker.
   Both checks would fail in a correctly configured Dfam 4.0 environment, producing a
   false-positive error and generating a useless `configure_dfam39.sh` script.

3. **`scripts/configure_dfam.sh`** — referenced Dfam 3.9 download URLs, `dfam39_full`
   filenames, and `perl ./configure`, all of which are obsolete.

### Files changed

**`rules/lib_construct.smk`** — `extract_repeatmasker_library` shell block:
- Replaced hardcoded path derivation with a two-tier detection:
  1. Search `$CONDA_PREFIX/share/famdb-*` for a standalone FamDB package (new, Dfam 4.0)
  2. Fall back to `share/RepeatMasker/Libraries/famdb/` (old, Dfam 3.9)
- Removed the `export PATH=...` line (famdb.py has been in `bin/` in all versions).

**`scripts/on_start_functions.py`** — EarlGrey configuration check block and helper:
- Replaced `library_path = rm_path.replace(...)` with the same two-tier `glob` detection.
- Extended h5 file check to match both `dfam40` and `dfam39_full` filenames.
- `.earlgrey.config.complete` marker is now a fast-path shortcut only (no longer required).
- Replaced `generate_dfam39_config_script()` with two new private helpers:
  - `_print_dfam_setup_instructions(library_path, famdb_shares)` — prints the correct
    setup instructions for either Dfam 4.0 (`download_dfam.py`) or Dfam 3.9 (curl steps).
  - `_generate_dfam39_config_script_legacy(library_path)` — generates the legacy shell
    script for Dfam 3.9 environments only.

**`scripts/configure_dfam.sh`** — rewritten as a reference guide:
- Describes the new `download_dfam.py` setup flow for Dfam 4.0.
- Preserves the old Dfam 3.9 manual steps as commented-out legacy notes.

**`conda/meta.yaml`** — version bumped to `0.2.0`; `earlgrey` minimum raised to `>=7.3.0`.

**`earlGreyParTEA`, `earlGreyParTEA_LibConstruct`, `earlGreyParTEA_AnnotationOnly`** —
`VERSION="0.2.0"`.

**`README.md`** — updated dependency badge, Configure EarlGrey setup section, Changes in
Latest Release section (v0.2.0 added; v0.1.9 moved to Previous Release).

### What did NOT change

- `repeatmasker_warmup` rule — the `CONS-*` dynamic discovery already matches `CONS-Dfam_4.0`.
- All RepeatMasker, RepeatModeler, and BuildDatabase invocation flags — unchanged.
- All Snakemake rule logic, clustering, annotation, saturation, shared/unique, BUSCO modules.
- Entry-point scripts and `generate_config.py` — contain no famdb-specific logic.

### Backward compatibility

Both old (≤7.2, Dfam 3.9) and new (7.3.0, Dfam 4.0) EarlGrey environments are supported.
The two-tier path detection falls back to the old embedded path when no `famdb-*` directory
exists in `share/`. The h5 file check matches either naming convention. The `.earlgrey.config.complete`
marker (if present from old setups) still triggers the fast-path skip.

### Verification checklist

- [X] With `earlgrey_730_test` active: `earlGreyParTEA_LibConstruct --dry-run` — confirm
      `[INFO] Found N Dfam partition file(s) in .../famdb-3.0.0/Libraries/famdb` in startup log

```
conda activate partea_20_test
earlGreyParTEA_LibConstruct --generate-config my_config.yaml --genome-dir /data/toby/testDIR/drosophila_testSet --output-dir /data/toby/testDIR/
earlGreyParTEA_LibConstruct -c my_config.yaml --threads 16 --dry-run
# this errored with: Please ensure EarlGrey (>=7.3.0) is properly installed and Dfam partitions have been downloaded via: download_dfam.py
download_dfam.py
# after running this, the dry run works
```

- [X] `extract_repeatmasker_library` step with `repeatmasker_species: "fungi"` — confirm
      `famdb.py -i .../famdb-3.0.0/Libraries/famdb/` in rule log
- [X] `repeatmasker_warmup` with `repeatmasker_species: "fungi"` — confirm
      `CONS-Dfam_4.0/fungi/` cache built successfully

```
# modify config to include RepeatMasker run with fungi
earlGreyParTEA_LibConstruct -c my_config.yaml --threads 16 --dry-run
earlGreyParTEA_LibConstruct -c my_config.yaml --threads 16
```

- [X] Full libconstruct run (4 genomes, no `repeatmasker_species`) — no regressions

