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