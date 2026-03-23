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